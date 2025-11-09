# lib_bpaf

Binary format for genomic sequence alignments with tracepoints.

## Features

- **Smart compression (recommended)**: Auto-selects optimal strategy based on data distribution
- **Multiple compression strategies**:
  - **Varint**: FastGA-aware delta + varint + zstd
  - **VarintRaw**: No delta encoding (test if delta hurts)
  - **Huffman**: FastGA-aware delta + adaptive Huffman + zstd
  - **HuffmanRaw**: No delta + adaptive Huffman (maximum entropy)
- **O(1) random access**: External index for instant record lookup
- **Tracepoint support**: Standard, Mixed, Variable, and FastGA representations
- **String deduplication**: Shared sequence name table

## Format

```
[Header] → [Codecs?] → [Records] → [StringTable]
```

### Header (6+ bytes)
- Magic: `BPAF` (4 bytes)
- Version: `1` (1 byte)
- Strategy: `0-5` (1 byte) - varint, varint-raw, huffman, huffman-raw, varint-huffman, varint-huffman-raw
- Record count: varint
- String count: varint

### Codecs (only for Huffman strategies 2-5)
- First value Huffman tree (symbol→bit mapping)
- Second value Huffman tree

### Records (per alignment)
- **PAF fields**: varints for positions, 1-byte strand/quality
- **Tracepoint metadata**: type (1 byte), complexity metric (1 byte), max_complexity (varint)
- **Tracepoints**: Compressed (first_values, second_values) pairs
- **Tags**: Optional key-value pairs

### String Table
- Deduplicated sequence names (length + UTF-8 bytes)

## Installation

```toml
[dependencies]
lib_bpaf = { git = "https://github.com/AndreaGuarracino/lib_bpaf" }
```

## Quick Start

### Read with random access

```rust
use lib_bpaf::BpafReader;

// Open with index for O(1) record access
let mut reader = BpafReader::open("alignments.bpaf")?;
println!("Total records: {}", reader.len());

// Jump to any record instantly
let record = reader.get_alignment_record(1000)?;
let (tracepoints, tp_type, _, _) = reader.get_tracepoints(1000)?;

match &tracepoints {
    TracepointData::Standard(tps) => println!("{} tracepoints", tps.len()),
    TracepointData::Fastga(tps) => println!("{} FastGA traces", tps.len()),
    _ => {}
}
```

### Read with file offsets (faster open)

```rust
// Skip index loading if you store offsets externally
let mut reader = BpafReader::open_without_index("alignments.bpaf")?;

// Access by byte offset
let offset = 123456;
let record = reader.get_alignment_record_at_offset(offset)?;
```

### Sequential iteration

```rust
for record in reader.iter_records() {
    let record = record?;
    println!("{} → {}", record.query_name_id, record.target_name_id);
}
```

### Compression

```rust
use lib_bpaf::{compress_paf, CompressionStrategy};

// Smart (recommended): Auto-select optimal strategy
// - Analyzes data distribution (top 3 symbol coverage)
// - Skewed data (>60%): Uses Huffman
// - Flat distribution: Uses Varint
compress_paf("alignments.paf", "alignments.bpaf", CompressionStrategy::Smart)?;

// Varint: FastGA-aware delta + varint + zstd
// - FastGA: raw num_differences (naturally small)
// - Standard: delta-encoded query_bases
compress_paf("alignments.paf", "alignments.bpaf", CompressionStrategy::Varint)?;

// Varint-Raw: No delta + varint + zstd
// - All types: raw first values (tests if delta hurts)
compress_paf("alignments.paf", "alignments.bpaf", CompressionStrategy::VarintRaw)?;

// Huffman: FastGA-aware delta + Huffman + zstd
// - FastGA: raw num_differences (exploits natural skew)
// - Standard: delta-encoded query_bases
compress_paf("alignments.paf", "alignments.bpaf", CompressionStrategy::Huffman)?;

// Huffman-Raw: No delta + Huffman + zstd
// - All types: raw first values (maximum entropy coding)
compress_paf("alignments.paf", "alignments.bpaf", CompressionStrategy::HuffmanRaw)?;
```

#### Huffman metrics (debug mode)

Tracepoint format varies by type:
- **Standard/Mixed/Variable**: `(query_bases, target_bases)` pairs
- **FastGA**: `(num_differences, target_bases)` pairs

Huffman compression encodes:
- **First values** (strategy-dependent):
  - `Huffman`: FastGA uses raw num_differences, Standard uses deltas
  - `HuffmanRaw`: All types use raw first values
- **Second values**: Target bases (all types, always raw)

Run with `RUST_LOG=debug` to see compression statistics:

```
[DEBUG] First value delta statistics:
  Unique symbols: 185
  Most frequent: 1 (45.2% of data)
  Top 3 symbols: 67.9% of data
  Code lengths: min=2, max=9, avg=4.56
  Theoretical bits/symbol: 8.00
  Huffman bits/symbol: 4.56
  Bit efficiency: 57.0%
[DEBUG] Second value statistics:
  Unique symbols: 180
  ...
```

**Metrics explained:**
- **Unique symbols**: Number of distinct values
- **Most frequent**: Most common value and its percentage
- **Top 3 symbols**: Coverage of three most common values
- **Theoretical bits**: Minimum bits with fixed-length encoding (log₂ symbols)
- **Huffman bits**: Weighted average bits per symbol
- **Bit efficiency**: Huffman/theoretical ratio (lower = better, >100% = worse than fixed)

**Strategy selection guide:**

**Recommended: Use `Smart`**
- Automatically analyzes your data and selects the best strategy
- No guesswork needed
- Logs the decision: `Smart: FastGA data with skew (top3=93.7%) → Huffman`

**Manual selection (if you want control):**

**For FastGA tracepoints:**
- High skew (top 3 >60%): Use `Huffman` (exploits natural num_differences distribution)
- Flat distribution: Use `Varint` (simpler, faster)

**For Standard tracepoints:**
- Start with `Varint` (FastGA-aware delta encoding)
- Try `VarintRaw` if delta encoding seems unhelpful (creates too many unique symbols)
- Use `Huffman`/`HuffmanRaw` only if you see high skew (top 3 >60%)

**General rules:**
- Skewed data (top 3 >60%, bit efficiency <80%): Huffman variants win by 10-20%
- Flat distributions: Varint variants win (simpler, faster, better compression)

### Index management

```rust
use lib_bpaf::{build_index, BpafIndex};

// Build index for random access
let index = build_index("alignments.bpaf")?;
index.save("alignments.bpaf.idx")?;

// Load existing index
let index = BpafIndex::load("alignments.bpaf.idx")?;
```

## Examples

```bash
# Build examples
cargo build --release --examples

# Show first 5 records
./target/release/examples/seek_demo alignments.bpaf

# O(1) random access demo
./target/release/examples/seek_demo alignments.bpaf 0 100 500 1000

# Offset-based access demo
./target/release/examples/offset_demo alignments.bpaf
```

## Index Format

`.bpaf.idx` file structure:
```
Magic:     BPAI (4 bytes)
Version:   1 (1 byte)
Count:     varint (number of records)
Offsets:   varint[] (byte positions)
```

Index enables O(1) random access without file scanning.

## License

MIT
