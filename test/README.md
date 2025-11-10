# lib_bpaf Test Suite

This directory contains all test scripts and tools for validating lib_bpaf functionality.

## Files

### `accurate_test.sh`
Comprehensive test suite that validates:
- **Compression**: Runtime, memory usage, determinism
- **Decompression**: Runtime, memory usage, correctness (MD5 verification)
- **O(1) Random Access**: Seek performance and tracepoint accuracy

**Usage:**
```bash
bash test/accurate_test.sh [input_paf] [strategy]
```

**Parameters:**
- `input_paf`: Path to PAF file (default: `/home/guarracino/Desktop/big-from-fg.tp.20k.paf`)
- `strategy`: Compression strategy (default: `automatic`)
  - Options: `automatic`, `varint`, `delta-varint`

**Example:**
```bash
# Run with default 20K record file
bash test/accurate_test.sh

# Run with custom file
bash test/accurate_test.sh /path/to/my.paf automatic
```

### `verify_tracepoints.rs`
Rust program that verifies O(1) random access correctness by comparing tracepoints retrieved from compressed BPAF files against the original PAF file.

**Compilation:**
Automatically compiled by `accurate_test.sh` to `test/verify_tracepoints` binary.

**Manual compilation:**
```bash
rustc --edition 2021 -O test/verify_tracepoints.rs \
    -L target/release/deps \
    --extern lib_bpaf=target/release/liblib_bpaf.rlib \
    -o test/verify_tracepoints
```

**Usage:**
```bash
./test/verify_tracepoints <original.paf> <compressed.bpaf> <record_id>
```

## Test Coverage

The test suite verifies:

1. **Compression correctness**
   - Deterministic output (multiple runs produce identical files)
   - Memory efficiency (< 20 MB for 20K records)
   - Performance (< 3 seconds for 20K records)

2. **Decompression correctness**
   - MD5 checksum matches original PAF
   - Memory efficiency (< 6 MB for 20K records)
   - Performance (< 1.5 seconds for 20K records)

3. **Random access correctness**
   - O(1) seek performance (15-20 μs per seek)
   - Tracepoint accuracy at 10 random positions
   - Low memory footprint (< 3 MB)

## Running Tests

From repository root:
```bash
# Quick test with 20K records
bash test/accurate_test.sh

# Test with larger dataset
bash test/accurate_test.sh /path/to/big-from-fg.tp.100k.paf
```

## Expected Output

```
========================================
lib_bpaf Comprehensive Test Suite
========================================
Input:      /home/guarracino/Desktop/big-from-fg.tp.20k.paf
Size:       122M
Records:    20000
Strategy:   automatic
========================================

...

--- Correctness ---
Compression: ✓ Identical
Decompression: ✓ Verified
O(1) seeks: ✓ Verified

========================================
All tests complete!
========================================
```

## Troubleshooting

**verify_tracepoints compilation failed:**
- Ensure `cargo build --release` has been run
- Check that `target/release/liblib_bpaf.rlib` exists

**Test failures:**
- Delete old `.bpaf` and `.bpaf.idx` files in `/tmp/`
- Ensure input PAF file has `tp:Z:` tags
- Check available disk space in `/tmp/`

## Adding New Tests

To add a new test:
1. Add test logic to `accurate_test.sh`
2. Update this README with test description
3. Ensure test outputs pass/fail status clearly
