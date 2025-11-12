#!/bin/bash
set -e

# Comprehensive test suite for lib_bpaf with proper tracepoint type testing
# Tests: encoding CIGAR→tracepoints, compression, decompression, O(1) seeks (Modes A/B)
# Measures: runtime, memory, correctness, compression ratios

REPO_DIR="/home/guarracino/Dropbox/git/lib_bpaf"
CIGZIP_DIR="/home/guarracino/Dropbox/git/cigzip"

# Source CIGAR-based PAF file
CIGAR_PAF="${1:-/home/guarracino/git/impg/hprcv2/data/hg002v1.1.pat.PanSN-vs-HG02818_mat_hprc_r2_v1.0.1.p95.Pinf.aln.paf.gz}"
MAX_COMPLEXITY="${2:-32}"
COMPLEXITY_METRIC="${3:-edit-distance}"
DISTANCE_METRIC="gap-affine"
PENALTIES="5,8,2"
NUM_RECORDS=20000
SEEK_ITERATIONS=1000
SEEK_PER_RECORD=100
PAF_SEEK_BIN="$REPO_DIR/target/release/examples/paf_seek_bench"

# Tracepoint types to test
TP_TYPES=("standard" "variable" "mixed")

# Verify input file exists
if [ ! -f "$CIGAR_PAF" ]; then
    echo "Error: CIGAR-based PAF file not found: $CIGAR_PAF"
    exit 1
fi

echo "========================================"
echo "lib_bpaf Comprehensive Test Suite"
echo "Proper Tracepoint Type Testing"
echo "========================================"
echo "Source:       $CIGAR_PAF"
echo "Records:      $NUM_RECORDS"
echo "Complexity:   $MAX_COMPLEXITY"
echo "Metric:       $COMPLEXITY_METRIC"
echo "Types:        ${TP_TYPES[*]}"
echo "========================================"
echo

# Arrays to store results
declare -A ENCODE_TIME COMPRESS_TIME COMPRESS_MEM COMPRESS_SIZE
declare -A DECOMP_TIME DECOMP_MEM
declare -A SEEK_TIME_A SEEK_TIME_B
declare -A VERIFIED

# Build cigzip
echo "=== Building cigzip ==="
if [ -d "$CIGZIP_DIR" ]; then
    cd "$CIGZIP_DIR"
    cargo build --release 2>&1 | tail -3
    CIGZIP="$CIGZIP_DIR/target/release/cigzip"

    if [ ! -f "$CIGZIP" ]; then
        echo "Error: cigzip binary not found after build"
        exit 1
    fi
else
    echo "Error: cigzip directory not found: $CIGZIP_DIR"
    exit 1
fi
echo "✓ cigzip built"
echo

# Build lib_bpaf test programs
echo "=== Building lib_bpaf test programs ==="
cd "$REPO_DIR"
cargo build --release 2>&1 | tail -3
cargo build --release --example paf_seek_bench 2>&1 | tail -3

# Build seek test program (Mode A - with index)
cat > /tmp/seek_test_mode_a.rs << 'RUST_EOF'
use std::env;
use std::time::Instant;
use lib_bpaf::BpafReader;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <bpaf> <record_id> <iterations>", args[0]);
        std::process::exit(1);
    }

    let bpaf_path = &args[1];
    let record_id: u64 = args[2].parse().expect("Invalid record_id");
    let iterations: usize = args[3].parse().expect("Invalid iterations");

    let mut reader = BpafReader::open(bpaf_path).expect("Failed to open BPAF");

    // Warmup
    for _ in 0..3 {
        let _ = reader.get_tracepoints(record_id).expect("Failed to fetch");
    }

    // Benchmark
    let mut total_seek_us = 0f64;
    for _ in 0..iterations {
        let seek_start = Instant::now();
        let _ = reader.get_tracepoints(record_id).expect("Failed to fetch");
        let seek_elapsed = seek_start.elapsed().as_micros() as f64;
        total_seek_us += seek_elapsed;
    }

    let avg_seek_us = total_seek_us / iterations as f64;
    println!("MODEA_SEEK {:.2}", avg_seek_us);
}
RUST_EOF

rustc --edition 2021 -O /tmp/seek_test_mode_a.rs \
    -L target/release/deps \
    --extern lib_bpaf=target/release/liblib_bpaf.rlib \
    -o /tmp/seek_test_mode_a 2>/dev/null

# Build seek test program (Mode B - standalone functions, fastest)
cat > /tmp/seek_test_mode_b.rs << 'RUST_EOF'
use std::env;
use std::time::Instant;
use std::fs::File;
use lib_bpaf::{BpafReader, read_standard_tracepoints_at_offset, 
               read_variable_tracepoints_at_offset, read_mixed_tracepoints_at_offset};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: {} <bpaf> <record_id> <iterations> <tp_type>", args[0]);
        std::process::exit(1);
    }

    let bpaf_path = &args[1];
    let record_id: u64 = args[2].parse().expect("Invalid record_id");
    let iterations: usize = args[3].parse().expect("Invalid iterations");
    let tp_type = &args[4];

    // Pre-compute offset and strategy using BpafReader (not timed)
    let mut reader = BpafReader::open(bpaf_path).expect("Failed to open BPAF");
    let offset = reader.get_tracepoint_offset(record_id).expect("Failed to get offset");
    let header = reader.header();
    let strategy = header.strategy().expect("Failed to get strategy");
    drop(reader);

    // Open file once for Mode B
    let mut file = File::open(bpaf_path).expect("Failed to open file");

    // Warmup
    for _ in 0..3 {
        match tp_type.as_str() {
            "standard" => { let _ = read_standard_tracepoints_at_offset(&mut file, offset, strategy); }
            "variable" => { let _ = read_variable_tracepoints_at_offset(&mut file, offset); }
            "mixed" => { let _ = read_mixed_tracepoints_at_offset(&mut file, offset); }
            _ => panic!("Invalid tp_type"),
        }
    }

    // Benchmark
    let mut total_seek_us = 0f64;
    for _ in 0..iterations {
        let seek_start = Instant::now();
        match tp_type.as_str() {
            "standard" => { let _ = read_standard_tracepoints_at_offset(&mut file, offset, strategy).expect("Failed"); }
            "variable" => { let _ = read_variable_tracepoints_at_offset(&mut file, offset).expect("Failed"); }
            "mixed" => { let _ = read_mixed_tracepoints_at_offset(&mut file, offset).expect("Failed"); }
            _ => panic!("Invalid tp_type"),
        }
        let seek_elapsed = seek_start.elapsed().as_micros() as f64;
        total_seek_us += seek_elapsed;
    }

    let avg_seek_us = total_seek_us / iterations as f64;
    println!("MODEB_SEEK {:.2}", avg_seek_us);
}
RUST_EOF

rustc --edition 2021 -O /tmp/seek_test_mode_b.rs \
    -L target/release/deps \
    --extern lib_bpaf=target/release/liblib_bpaf.rlib \
    -o /tmp/seek_test_mode_b 2>/dev/null

echo "✓ Build complete"
echo

# Extract sample records
echo "=== Extracting $NUM_RECORDS records from CIGAR PAF ==="
if [[ "$CIGAR_PAF" == *.gz ]]; then
    zcat "$CIGAR_PAF" | head -n "$NUM_RECORDS" > /tmp/test.cigar.paf
else
    head -n "$NUM_RECORDS" "$CIGAR_PAF" > /tmp/test.cigar.paf
fi
EXTRACTED=$(wc -l < /tmp/test.cigar.paf)
SIZE=$(du -h /tmp/test.cigar.paf | cut -f1)
SIZE_BYTES=$(stat -f%z /tmp/test.cigar.paf 2>/dev/null || stat -c%s /tmp/test.cigar.paf)
echo "Extracted: $EXTRACTED records"
echo "Size:      $SIZE ($SIZE_BYTES bytes)"
echo
echo

echo "=== Benchmarking PAF BGZF CIGAR seeks ==="
PAF_SEEK_TIME="N/A"
PAF_INDEXED=0
if [ -x "$PAF_SEEK_BIN" ]; then
    PAF_SEEK_OUTPUT=$("$PAF_SEEK_BIN" "$CIGAR_PAF" "$EXTRACTED" 10 "$SEEK_PER_RECORD" 2>&1 || true)
    echo "$PAF_SEEK_OUTPUT"
    if [[ $PAF_SEEK_OUTPUT =~ PAF_INDEXED[[:space:]]+([0-9]+) ]]; then
        PAF_INDEXED="${BASH_REMATCH[1]}"
    fi
    if [[ $PAF_SEEK_OUTPUT =~ PAF_SEEK[[:space:]]+([0-9.]+) ]]; then
        PAF_SEEK_TIME="${BASH_REMATCH[1]}"
    fi
else
    echo "Warning: paf_seek_bench binary not found at $PAF_SEEK_BIN"
fi
echo
echo

# Test each tracepoint type
for TP_TYPE in "${TP_TYPES[@]}"; do
    echo "###################################################################"
    echo "# TESTING TRACEPOINT TYPE: ${TP_TYPE^^}"
    echo "###################################################################"
    echo

    # TEST 1: Encode CIGAR → Tracepoints
    echo "========================================"
    echo "TEST 1 ($TP_TYPE): Encode CIGAR → Tracepoints"
    echo "========================================"
    
    START_TIME=$(date +%s.%N)
    $CIGZIP encode --paf /tmp/test.cigar.paf --threads 1 --type "$TP_TYPE" \
        --max-complexity "$MAX_COMPLEXITY" --complexity-metric "$COMPLEXITY_METRIC" \
        > /tmp/test.$TP_TYPE.tp.paf 2>/dev/null
    END_TIME=$(date +%s.%N)
    
    ENCODE_TIME[$TP_TYPE]=$(echo "$END_TIME - $START_TIME" | bc -l)
    TP_SIZE=$(du -h /tmp/test.$TP_TYPE.tp.paf | cut -f1)
    TP_SIZE_BYTES=$(stat -f%z /tmp/test.$TP_TYPE.tp.paf 2>/dev/null || stat -c%s /tmp/test.$TP_TYPE.tp.paf)
    
    echo "Encode time: ${ENCODE_TIME[$TP_TYPE]}s"
    echo "Output size: $TP_SIZE ($TP_SIZE_BYTES bytes)"
    echo

    # TEST 2: Compress Tracepoints → BPAF
    echo "========================================"
    echo "TEST 2 ($TP_TYPE): Compress Tracepoints → BPAF"
    echo "========================================"
    
    START_TIME=$(date +%s.%N)
    rm -f /tmp/test.$TP_TYPE.bpaf.idx
    COMPRESS_LOG=$(env RUST_LOG=info $CIGZIP compress -i /tmp/test.$TP_TYPE.tp.paf \
        -o /tmp/test.$TP_TYPE.bpaf \
        --type "$TP_TYPE" \
        --max-complexity "$MAX_COMPLEXITY" \
        --complexity-metric "$COMPLEXITY_METRIC" \
        --distance "$DISTANCE_METRIC" \
        --penalties "$PENALTIES" 2>&1)
    if ! echo "$COMPRESS_LOG" | grep -E "INFO|Compressed"; then
        echo "$COMPRESS_LOG"
    fi
    END_TIME=$(date +%s.%N)
    
    BPAF_SIZE=$(du -h /tmp/test.$TP_TYPE.bpaf | cut -f1)
    BPAF_SIZE_BYTES=$(stat -f%z /tmp/test.$TP_TYPE.bpaf 2>/dev/null || stat -c%s /tmp/test.$TP_TYPE.bpaf)
    COMPRESS_SIZE[$TP_TYPE]=$BPAF_SIZE_BYTES
    
    echo "Compress time: N/As"
    echo "Output size:   $BPAF_SIZE ($BPAF_SIZE_BYTES bytes)"
    echo

    # TEST 3: Decompress BPAF → Tracepoints
    echo "========================================"
    echo "TEST 3 ($TP_TYPE): Decompress BPAF → Tracepoints"
    echo "========================================"
    
    START_TIME=$(date +%s.%N)
    $CIGZIP decompress -i /tmp/test.$TP_TYPE.bpaf \
        -o /tmp/test.$TP_TYPE.decomp.paf
    END_TIME=$(date +%s.%N)
    
    DECOMP_TIME[$TP_TYPE]=$(echo "$END_TIME - $START_TIME" | bc -l)
    DECOMP_MEM[$TP_TYPE]="N/A"
    
    echo "Decomp time:  ${DECOMP_TIME[$TP_TYPE]}s"
    echo "Memory:       ${DECOMP_MEM[$TP_TYPE]} KB"

    # Normalize floats to 3 decimal places and compute MD5
    perl -pe 's/(\d+\.\d{3})\d*/$1/g' /tmp/test.$TP_TYPE.tp.paf | md5sum > /tmp/orig.md5
    perl -pe 's/(\d+\.\d{3})\d*/$1/g' /tmp/test.$TP_TYPE.decomp.paf | md5sum > /tmp/decomp.md5

    ORIG_MD5=$(cut -d' ' -f1 /tmp/orig.md5)
    DECOMP_MD5=$(cut -d' ' -f1 /tmp/decomp.md5)

    if [ "$ORIG_MD5" = "$DECOMP_MD5" ]; then
        echo "✓ Tracepoints verified: Perfect match (MD5: $ORIG_MD5, normalized to 3 decimal places for all floats)"
        VERIFIED[$TP_TYPE]="✓ Perfect match"
    else
        echo "✗ Tracepoints verification FAILED"
        echo "  Original MD5:      $ORIG_MD5"
        echo "  Decompressed MD5:  $DECOMP_MD5"
        VERIFIED[$TP_TYPE]="✗ MISMATCH"
    fi
    echo

    # TEST 4: O(1) Random Access
    echo "========================================"
    echo "TEST 4 ($TP_TYPE): O(1) Random Access"
    echo "========================================"
    echo "Testing $SEEK_ITERATIONS seeks at different positions..."
    echo

    # Mode A: BpafReader with index
    TOTAL_TIME_A=0
    COUNT_A=0
    for i in $(seq 0 10 $((EXTRACTED - 1))); do
        OUTPUT=$(/tmp/seek_test_mode_a /tmp/test.$TP_TYPE.bpaf "$i" "$SEEK_PER_RECORD" 2>/dev/null)
        if [[ $OUTPUT =~ MODEA_SEEK[[:space:]]+([0-9.]+) ]]; then
            TIME="${BASH_REMATCH[1]}"
            TOTAL_TIME_A=$(echo "$TOTAL_TIME_A + $TIME" | bc -l)
            COUNT_A=$((COUNT_A + 1))
        fi
    done

    if [ "$COUNT_A" -gt 0 ]; then
        SEEK_TIME_A[$TP_TYPE]=$(echo "scale=2; $TOTAL_TIME_A / $COUNT_A" | bc -l)
    else
        SEEK_TIME_A[$TP_TYPE]="N/A"
    fi

    # Mode B: Standalone functions (fastest)
    TOTAL_TIME_B=0
    COUNT_B=0
    for i in $(seq 0 10 $((EXTRACTED - 1))); do
        OUTPUT=$(/tmp/seek_test_mode_b /tmp/test.$TP_TYPE.bpaf "$i" "$SEEK_PER_RECORD" "$TP_TYPE" 2>/dev/null)
        if [[ $OUTPUT =~ MODEB_SEEK[[:space:]]+([0-9.]+) ]]; then
            TIME="${BASH_REMATCH[1]}"
            TOTAL_TIME_B=$(echo "$TOTAL_TIME_B + $TIME" | bc -l)
            COUNT_B=$((COUNT_B + 1))
        fi
    done

    if [ "$COUNT_B" -gt 0 ]; then
        SEEK_TIME_B[$TP_TYPE]=$(echo "scale=2; $TOTAL_TIME_B / $COUNT_B" | bc -l)
    else
        SEEK_TIME_B[$TP_TYPE]="N/A"
    fi

    echo "Average seek time (Mode A - BpafReader with index):     ${SEEK_TIME_A[$TP_TYPE]} μs"
    echo "Average seek time (Mode B - Standalone functions):      ${SEEK_TIME_B[$TP_TYPE]} μs"
    echo
    echo
done

# Final summary
echo "###################################################################"
echo "# FINAL SUMMARY - ALL TRACEPOINT TYPES"
echo "###################################################################"
echo
echo "Source:  $CIGAR_PAF"
echo "Sample:  $EXTRACTED records, $SIZE ($SIZE_BYTES bytes)"
echo

# Summary table
printf "╔═══════════╦══════════════╦═══════════════╦════════════════╦═══════════════╦═════════════╦═════════════╗\n"
printf "║ %-9s ║ %-12s ║ %-13s ║ %-14s ║ %-13s ║ %-11s ║ %-11s ║\n" \
    "Type" "Encode (s)" "Compress (s)" "Decomp (s)" "Size (bytes)" "Seek A (μs)" "Seek B (μs)"
printf "╠═══════════╬══════════════╬═══════════════╬════════════════╬═══════════════╬═════════════╬═════════════╣\n"

for TP_TYPE in "${TP_TYPES[@]}"; do
    printf "║ %-9s ║ %12.3f ║ %13s ║ %14s ║ %13s ║ %11s ║ %11s ║\n" \
        "$TP_TYPE" \
        "${ENCODE_TIME[$TP_TYPE]}" \
        "N/A" \
        "${DECOMP_TIME[$TP_TYPE]}" \
        "${COMPRESS_SIZE[$TP_TYPE]}" \
        "${SEEK_TIME_A[$TP_TYPE]}" \
        "${SEEK_TIME_B[$TP_TYPE]}"
done

printf "╚═══════════╩══════════════╩═══════════════╩════════════════╩═══════════════╩═════════════╩═════════════╝\n"
echo

echo "Data Integrity:"
for TP_TYPE in "${TP_TYPES[@]}"; do
    printf "  %-10s: %s\n" "$TP_TYPE" "${VERIFIED[$TP_TYPE]}"
done
echo

echo "Compression Ratios (CIGAR PAF → BPAF):"
for TP_TYPE in "${TP_TYPES[@]}"; do
    RATIO=$(echo "scale=2; $SIZE_BYTES / ${COMPRESS_SIZE[$TP_TYPE]}" | bc -l)
    printf "  %-10s: %sx (CIGAR: %s bytes → BPAF: %s bytes)\n" \
        "$TP_TYPE" "$RATIO" "$SIZE_BYTES" "${COMPRESS_SIZE[$TP_TYPE]}"
done
echo

echo "PAF BGZF (cg:Z:) seek benchmark:"
echo "  Records indexed: $PAF_INDEXED"
echo "  Average seek time: ${PAF_SEEK_TIME} μs (stride 10, ${SEEK_PER_RECORD} iterations)"
echo

echo "Seek Mode Legend:"
echo "  Mode A: BpafReader::open() + get_tracepoints()              (with index, general use)"
echo "  Mode B: read_*_tracepoints_at_offset(File)  (standalone functions, ultimate performance)"
echo

echo "========================================"
echo "All tests complete!"
echo "========================================"
