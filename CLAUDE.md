# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`tpa` is a Rust library for binary compression and random access of genomic sequence alignments with tracepoints. It provides O(1) random record access through external indexing and supports multiple compression strategies optimized for genomic data.

## Build & Test Commands

### Building
```bash
# Build library
cargo build --release

# Build specific binary
cargo build --release --bin tpa-analyze

# Build all examples
cargo build --release --examples
```

### Testing
```bash
# Run Rust unit/integration tests
cargo test

# Comprehensive test suite (tests dual-strategy combinations across multiple datasets)
# Usage: ./scripts/run_all_tests.sh [num_records] [output_base] [paf1] [paf2] ... [pafN]
./scripts/run_all_tests.sh                    # Use all defaults (50 records, 3 default PAF files)
./scripts/run_all_tests.sh 100                # Test 100 records per file
./scripts/run_all_tests.sh 200 /tmp/results   # Custom output directory
./scripts/run_all_tests.sh 50 /tmp/out file1.paf.gz file2.paf  # Custom files

# Test single file with auto-detection (CIGAR or tracepoint PAF)
./scripts/comprehensive_test.sh input.paf.gz [output_dir] [max_complexity] [metric] [records]

# Example: Test with 1000 records
./scripts/comprehensive_test.sh data.paf.gz /tmp/output 32 edit-distance 1000

# Check if automatic mode selected optimal strategies
./scripts/check-automatic.py [all_results.tsv]
```

The test suite produces:
- Individual test reports: `output_dir/{dataset}/test_report.md`
- Aggregated results: `output_dir/FINAL_REPORT.md`
- TSV data: `output_dir/all_results.tsv` (32 columns including dual strategies, per-stream layers, compression ratios, seek performance, verification status)

### Running Examples
```bash
# Analyze TPA file structure
./target/release/tpa-analyze file.tpa

# View TPA contents as PAF (decompress to stdout)
./target/release/tpa-view file.tpa

# Show which strategies/layers were selected (useful for automatic mode)
./target/release/tpa-view --strategies file.tpa

# Seek performance demo
./target/release/examples/seek_demo file.tpa [record_ids...]

# Offset-based access demo
./target/release/examples/offset_demo file.tpa
```

### Benchmark Tool

Unified seek performance benchmark located in `scripts/`:

```bash
# Build unified benchmark (auto-detects format)
rustc scripts/seek_bench.rs --edition 2021 \
    --extern tpa=target/release/libtpa.rlib \
    --extern noodles=target/release/deps/libnoodles*.rlib \
    -L target/release/deps -o /tmp/seek_bench

# Benchmark any format - auto-detects CIGAR PAF, tracepoint PAF, or TPA
/tmp/seek_bench input.paf.gz 10000 100 100
/tmp/seek_bench input.tpa 10000 100 100
```

## Architecture

### Core Components

**src/lib.rs** - Public API surface
- `TpaReader`: Main reader with O(1) random access via index
- `compress_paf_to_tpa()`: Compression entry point (uses `CompressionConfig` builder)
- `read_*_tracepoints_at_offset()`: Fast standalone seek functions (no TpaReader overhead)
- `build_index()`: Creates `.tpa.idx` files for random access

**src/format.rs** - Data structures and serialization
- `CompressionStrategy` enum: 19 concrete strategies + 1 meta-strategy (`Automatic(level, sample_size)`)
- `CompressionLayer` enum: Zstd, Bgzip, Nocomp - passed explicitly through API
- `TpaHeader`: File metadata (magic, version=1, first/second compression layers, first/second strategy codes, counts, tracepoint type, distance params)
- `TpaFooter`: Crash-safety footer with record/string counts for validation
- `AlignmentRecord`: PAF fields (positions, strand, quality, name IDs)
- `TracepointData` enum: Standard/Variable/Mixed/FastGA representations
- **Dual strategies**: Header stores separate codes for 1st and 2nd values (first_strategy_code, second_strategy_code)

**src/binary.rs** - Binary encoding/decoding
- `SmartDualAnalyzer`: Tests every concrete strategy × layer per stream (19×3 each) and picks independent winners for first/second values
- Strategy-specific encode/decode for tracepoints (19 strategies)
- Varint encoding with delta/zigzag transforms
- Compression layer handling (Zstd/Bgzip/Nocomp)
- Rice and Huffman entropy coding

**src/hybrids.rs** - Advanced compression strategies
- FastPFOR: Patched Frame-of-Reference with exceptions
- Cascaded: Dictionary → RLE chains for low-cardinality data
- Simple8bFull: 16-selector word packing
- SelectiveRLE: Run-length preprocessing with bitmap positions

**src/reader.rs** - High-level reader API
- `TpaReader`: Main reader struct with header, index, string table
- `read_*_tracepoints_at_offset()`: Fast standalone seek functions
- `RecordIterator`: Sequential record iteration

**src/index.rs** - Index building and management
- `TpaIndex`: Index struct with record byte offsets
- `build_index()`: Scans TPA file to create index
- Index save/load for `.tpa.idx` files

**src/utils.rs** - Varint utilities and file helpers
- `write_varint()`, `read_varint()`: Basic varint I/O
- `encode_zigzag()`, `decode_zigzag()`: Signed-to-unsigned transforms
- `open_paf_reader()`: Opens PAF files (handles .gz compression)

### Binary Format

```
[Header] → [StringTable] → [Records...] → [Footer]
```

**Common Prefix (shared by Header and Footer):**
- Magic: "TPA\0" (4 bytes, null-terminated for C compatibility and 32-bit alignment)
- Version: 1 (1 byte)
- Record count: varint
- String count: varint

**Header (Version 1):** Common prefix + header-specific fields:
- First strategy+layer: packed byte (6 bits strategy code 0-18, 2 bits layer)
- Second strategy+layer: packed byte (6 bits strategy code 0-18, 2 bits layer)
- Tracepoint type: 1 byte (Standard/Variable/Mixed/FastGA)
- Complexity metric: 1 byte
- Max complexity/spacing: varint
- Distance parameters: serialized for lib_wfa2

**Footer:** Common prefix + footer_length:
- Footer length: 4 bytes (little-endian u32)

**Index Format (.tpa.idx):**
- Magic: "TPAI" (4 bytes)
- Version: 1 (1 byte)
- Record count: varint
- Byte offsets: varint array (one per record)

## Compression Strategy Selection

### Automatic Strategy
`Automatic(level, sample_size)` - Tests every concrete strategy × layer per stream (19×3=57 combos each), selects independent winners for first/second values.
- `sample_size=1000` (default): Samples first 1,000 records for fast analysis
- `sample_size=0`: Analyzes entire file for best compression (slower)
- `sample_size=N`: Custom sample size

Dual encoding (different strategies for 1st/2nd values) is achieved via `CompressionConfig::dual_strategy()`.

### Concrete Strategies (19 total)
All 19 strategies are considered during automatic analysis:

**High performers:**
1. **2d-delta**: Best for CIGAR-derived alignments (exploits query/target correlation)
2. **stream-vbyte**: SIMD-friendly byte-aligned encoding
3. **simple8b-full**: Complete Simple8b with 16 packing modes + RLE
4. **zigzag-delta**: General-purpose fallback
5. **cascaded**: Multi-level encoding for low-cardinality data

**Specialized:**
- **raw**: Low-complexity data
- **selective-rle**: High-repetition blocks
- **dictionary**: Repeated delta patterns
- **rice**, **huffman**: Entropy coding alternatives

## Key Implementation Details

### Tracepoint Encoding
All strategies preserve **byte-aligned varint encoding** for O(1) random access:
- First values and second values stored separately
- Each value type can use a different compression strategy (dual encoding)
- Delta encoding (when used) maintains position-independence per record
- Compression layers (Zstd/Bgzip/Nocomp) stored per stream and applied after varint encoding
- Layer parameter passed explicitly through all API functions (no thread-local state)

### Random Access Pattern
1. Load index (.tpa.idx) to get byte offsets
2. Seek to record offset in .tpa file
3. Decode tracepoint data using stored strategies from header (first_strategy_code, second_strategy_code)
4. No need to decode previous records

**Dual strategies**: First and second values decoded independently with their respective strategies

### Float Normalization (Testing)
Test suite uses `scripts/normalize_paf.py` for verification:
- Truncates floats to 3 decimals: `0.993724` → `0.993`
- Handles edge cases: `.0549` → `0.054`, `0` → `0.000`

## Dependencies

**Core:**
- `zstd`: Compression layer
- `tracepoints`: Tracepoint representations and CIGAR→tracepoint conversion
- `lib_wfa2`: WFA alignment distance parameters

**Testing:**
- `noodles`: BGZF decompression for .paf.gz files
- Python 3 (scripts/normalize_paf.py), md5sum, /usr/bin/time
- **Rust toolchain**: cargo and rustc must be in PATH (sourced via ~/.cargo/env)

## Troubleshooting

### Test script dies at "Building seek test programs"

**Symptom**: `comprehensive_test.sh` hangs or fails silently at "=== Building seek test programs ==="

**Root cause**: `rustc` command not found in PATH

**Solution**: The script now automatically sources `~/.cargo/env` to set up Rust environment. If this fails:
1. Verify Rust is installed: `which rustc`
2. Check that `~/.cargo/env` exists and is valid
3. Manually source before running tests: `source ~/.cargo/env && ./scripts/run_all_tests.sh`

**Prevention**: The script includes proper error checking for rustc compilation failures (scripts/comprehensive_test.sh)

## Recent Development (2025-11-28)

### Directory Restructure and Cleanup

**Changes Made:**
- Moved test/utility scripts from `test/` to `scripts/` directory
- Removed `tpa-viz` binary (visualization feature removed)
- Added `scripts/check-automatic.py` to verify automatic mode selects optimal strategies

### Header/Footer Format Alignment

**Changes Made:**
- Reordered header fields to match footer: `[magic][version][num_records][num_strings]` prefix is now identical
- Added `write_common_prefix()` and `read_common_prefix()` helper functions in `src/format.rs`
- Both `TpaHeader` and `TpaFooter` now use shared helpers for the common prefix
- Test scripts updated: Added `AUTO_STRATEGIES` array to `scripts/comprehensive_test.sh` to fix automatic mode testing in dual mode

**Benefits:**
- ~25 lines of code reduced through deduplication
- Consistent validation logic for magic/version/counts
- Cleaner separation between common metadata and format-specific fields
- Easier header/footer cross-validation

### Repository Rename: lib_bpaf → tpa

**Changes Made:**
- GitHub repository renamed from `lib_bpaf` to `tpa`
- Updated test scripts to use `tpa_` prefixes instead of `bpaf_`
- TSV column names updated: `bpaf_size_bytes` → `tpa_size_bytes`, `ratio_tp_to_bpaf` → `ratio_tp_to_tpa`, etc.
- Default output directory: `scripts/tpa_all_tests/`
- Rust API already used `TpaReader` (no code changes needed)
- Renamed internal structs: `BinaryPafHeader` → `TpaHeader`, `BinaryPafFooter` → `TpaFooter`

---

## Earlier Development (2025-11-26)

### Code Consolidation and API Cleanup

**Changes Made:**

1. **New CompressionConfig API**: Replaced 4 wrapper functions with single `compress_paf_to_tpa()` using builder pattern
   ```rust
   compress_paf_to_tpa("input.paf", "output.tpa",
       CompressionConfig::new()
           .strategy(CompressionStrategy::ZigzagDelta(3))
           .layer(CompressionLayer::Zstd)
           .from_cigar()
   )?;
   ```

2. **Consolidated zigzag encoding**: Moved `encode_zigzag()` and `decode_zigzag()` to utils.rs for sharing across modules

3. **Consolidated decode match arms**: Combined ZigzagDelta, TwoDimDelta, and OffsetJoint into single match arm (they share identical decode logic)

4. **Removed dead code**: Deleted unused `is_strategy_name()` function, removed incorrect `#[allow(dead_code)]` attributes

**Lines reduced**: 80 lines (5155 → 5075)

---

## Development (2025-11-17)

### Dual Compression Strategy Implementation

**Major Features Completed:**

1. **Dual Strategy Support**: Separate compression strategies for 1st and 2nd values in tracepoint pairs
   - Header stores `first_strategy_code`, `second_strategy_code`, and two compression-layer bytes (version 1 format)
   - Compression layers recorded independently (first/second streams can mix Zstd/Bgzip/Nocomp)
   - When codes differ: creates Dual strategy; when same: single strategy

2. **Unified Automatic Mode**: `Automatic(level, sample_size)` tests 57 combinations per stream (19 strategies × 3 layers)
   - Returns optimal tuple: (first_strategy, first_layer, second_strategy, second_layer)
   - sample_size=1000 (default) for fast sampling; sample_size=0 for entire file
   - Includes Rice and Huffman along with all other concrete codecs

3. **Thread-Local State Removal**: Replaced with explicit parameter passing
   - `layer` parameter added to all public APIs
   - Thread-safe: multiple threads can compress with different layers simultaneously
   - Cleaner architecture: explicit data flow throughout call stack

4. **Tool Consolidation**:
   - Unified seek benchmark: `scripts/seek_bench.rs` auto-detects format (CIGAR PAF, tracepoint PAF, or TPA)
   - Removed redundant tools: compress_paf, tpa_header, separate benchmark tools

5. **Test Suite Enhancements**:
   - Tests all 3,250 combinations per file (3,249 explicit dual combos + 1 automatic)
   - Dual mode: `comprehensive_test.sh` with TEST_MODE=dual parameter
   - Uses `tpa-view --strategies` to extract and track which strategies and layers automatic mode selected
   - Plot visualization: 4×1 vertical layout with improved readability

**Key Files Modified:**
- src/format.rs: Dual strategy support, CompressionLayer enum (21 total strategy variants: 2 meta + 19 concrete)
- src/binary.rs: SmartDualAnalyzer with per-stream 57-combo testing (19 strategies × 3 layers)
- src/lib.rs: Dual strategy APIs, explicit layer parameter
- scripts/comprehensive_test.sh: Dual mode support, strategy tracking
- scripts/run_all_tests.sh: Flexible parameterization

**Dataset-Specific Insights:**
- CIGAR-derived data (p95, sweepga): TwoDimDelta often optimal for 2nd values
- Native tracepoint data (bigfg): Simple strategies (Raw, Simple8) often win
- Compression layer choice is dataset-dependent (Zstd vs Nocomp varies)

## Common Patterns

### Adding a New Compression Strategy

1. Add enum variant to `CompressionStrategy` in `src/format.rs`:
   ```rust
   pub enum CompressionStrategy {
       // ...
       MyNewStrategy(i32), // parameter is zstd level
   }
   ```

2. Add encode/decode in `src/binary.rs`:
   ```rust
   CompressionStrategy::MyNewStrategy(_) => {
       // Encode logic - must produce byte-aligned varint stream
   }
   ```

3. Add to strategy list in `src/binary.rs` (~line 160)

4. Add to test script `scripts/comprehensive_test.sh` STRATEGIES array

5. Update `src/lib.rs` strategy resolution if meta-strategy

### Reading Performance-Critical Code

When tracepoints must be decoded in tight loops:
- Use `read_*_tracepoints_at_offset()` functions (src/lib.rs)
- Pre-compute offsets from index
- Avoid `TpaReader` overhead
- See `examples/offset_demo.rs` for pattern

### Verification Failures

If round-trip tests fail:
1. Check float normalization: `python3 ./scripts/normalize_paf.py file.paf | head`
2. Verify strategy preserves byte-alignment
3. Test with smaller dataset to isolate issue
4. Check that delta encoding is position-independent per record
