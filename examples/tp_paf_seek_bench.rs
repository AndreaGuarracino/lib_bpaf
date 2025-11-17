/// Benchmark tool for testing O(1) random access performance on bgzipped tracepoint PAF files
///
/// This tool:
/// 1. Builds an index of tracepoint positions (using BGZF virtual positions)
/// 2. Performs random seeks to test access performance
/// 3. Reports statistics on seek times (average, stddev, min, max)
///
/// Usage: tp_paf_seek_bench <tp.paf.bgz> <num_records> <num_positions> <iterations_per_position>
///
/// Example: tp_paf_seek_bench alignments.tp.paf.gz 10000 100 100
///   - Index first 10000 records
///   - Test 100 random positions
///   - Do 100 iterations per position

use noodles::bgzf::io::Reader;
use noodles::bgzf::VirtualPosition;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Read, Write};
use std::path::Path;
use std::time::Instant;

#[derive(Clone, Copy, Debug)]
struct TracepointEntry {
    record_id: usize,
    pos: VirtualPosition,
    len: usize,
}

impl TracepointEntry {
    fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.record_id.to_le_bytes());
        bytes.extend_from_slice(&u64::from(self.pos).to_le_bytes());
        bytes.extend_from_slice(&self.len.to_le_bytes());
        bytes
    }

    fn from_bytes(bytes: &[u8]) -> Result<Self, Box<dyn std::error::Error>> {
        if bytes.len() < 24 {
            return Err("Invalid index entry size".into());
        }
        let record_id = usize::from_le_bytes(bytes[0..8].try_into()?);
        let pos_u64 = u64::from_le_bytes(bytes[8..16].try_into()?);
        let pos = VirtualPosition::from(pos_u64);
        let len = usize::from_le_bytes(bytes[16..24].try_into()?);
        Ok(Self { record_id, pos, len })
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!(
            "Usage: {} <tp.paf.bgz> <num_records> <num_positions> <iterations>",
            args[0]
        );
        eprintln!();
        eprintln!("Arguments:");
        eprintln!("  <tp.paf.bgz>        : Path to bgzipped tracepoint PAF file");
        eprintln!("  <num_records>       : Number of records to index (0 = all)");
        eprintln!("  <num_positions>     : Number of random positions to test");
        eprintln!("  <iterations>        : Iterations per position");
        eprintln!();
        eprintln!("Example:");
        eprintln!("  {} alignments.tp.paf.gz 10000 100 100", args[0]);
        std::process::exit(1);
    }

    let paf_path = &args[1];
    let num_records: usize = args[2]
        .parse()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid num_records"))?;
    let num_positions: usize = args[3]
        .parse()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid num_positions"))?;
    let iterations: usize = args[4]
        .parse()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid iterations"))?;

    if num_positions == 0 || iterations == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "num_positions and iterations must be greater than zero",
        )
        .into());
    }

    println!("=== Tracepoint PAF Seek Benchmark ===");
    println!("File: {}", paf_path);

    let index_path = format!("{}.idx", paf_path);
    let entries = if Path::new(&index_path).exists() {
        println!("Loading existing index...");
        let start = Instant::now();
        let mut entries = load_index(&index_path)?;
        let load_time = start.elapsed();

        // Truncate if num_records is specified and less than indexed
        if num_records > 0 && num_records < entries.len() {
            entries.truncate(num_records);
        }

        println!(
            "Loaded {} records in {:.2}s",
            entries.len(),
            load_time.as_secs_f64()
        );
        entries
    } else {
        println!("Building index...");
        let start = Instant::now();
        let entries = build_tracepoint_index(paf_path, num_records)?;
        let index_time = start.elapsed();

        if entries.is_empty() {
            println!("No tracepoint records found!");
            return Ok(());
        }

        println!(
            "Indexed {} records in {:.2}s",
            entries.len(),
            index_time.as_secs_f64()
        );

        // Save index for future use
        println!("Saving index to {}...", index_path);
        save_index(&entries, &index_path)?;

        entries
    };

    if entries.is_empty() {
        println!("No tracepoint records found!");
        return Ok(());
    }

    let test_limit = if num_records == 0 || num_records > entries.len() {
        entries.len()
    } else {
        num_records
    };

    println!();

    // Select random positions to test
    let test_positions = select_random_positions(&entries, num_positions);
    println!(
        "Testing {} random positions × {} iterations each",
        test_positions.len(),
        iterations
    );
    println!();

    let stats = benchmark_tracepoint_seeks(paf_path, &entries, &test_positions, iterations)?;

    println!("=== Results ===");
    println!("Records available:   {}", entries.len());
    println!("Records tested:      {}", test_limit);
    println!("Positions tested:    {}", test_positions.len());
    println!("Iterations/position: {}", iterations);
    println!("Total seeks:         {}", test_positions.len() * iterations);
    println!();
    println!("Seek Time (microseconds):");
    println!("  Average:  {:.2} μs", stats.avg);
    println!("  Std Dev:  {:.2} μs", stats.stddev);
    println!("  Median:   {:.2} μs", stats.median);
    println!("  Min:      {:.2} μs", stats.min);
    println!("  Max:      {:.2} μs", stats.max);
    println!("  P95:      {:.2} μs", stats.p95);
    println!("  P99:      {:.2} μs", stats.p99);
    println!();
    println!("Success rate: {:.2}%", stats.success_rate * 100.0);
    println!();

    // Check if O(1) is achieved (< 1ms for in-memory access)
    if stats.avg < 1000.0 {
        println!("✓ O(1) performance achieved (mean < 1ms)");
    } else {
        println!("⚠ Slow seeks detected (mean >= 1ms)");
    }

    // Check consistency (low variance means predictable O(1))
    let cv = (stats.stddev / stats.avg) * 100.0; // Coefficient of variation
    println!("  Coefficient of variation: {:.1}%", cv);

    if cv < 20.0 {
        println!("✓ Consistent O(1) performance (CV < 20%)");
    } else if cv < 50.0 {
        println!("⚠ Moderate variance in seek times (CV < 50%)");
    } else {
        println!("⚠ High variance in seek times (CV >= 50%)");
    }

    // Output in machine-readable format for test scripts
    println!();
    println!("# Machine-readable output:");
    println!("TP_PAF_INDEXED {}", entries.len());
    println!("TP_PAF_SEEK_AVG {:.2}", stats.avg);
    println!("TP_PAF_SEEK_STDDEV {:.2}", stats.stddev);
    println!("TP_PAF_SUCCESS_RATE {:.4}", stats.success_rate);

    Ok(())
}

#[derive(Debug)]
struct SeekStats {
    avg: f64,
    stddev: f64,
    median: f64,
    min: f64,
    max: f64,
    p95: f64,
    p99: f64,
    success_rate: f64,
}

fn build_tracepoint_index(
    paf_path: &str,
    limit: usize,
) -> Result<Vec<TracepointEntry>, Box<dyn std::error::Error>> {
    let file = File::open(paf_path)?;
    let mut reader = Reader::new(file);
    let mut line = Vec::new();
    let mut entries = Vec::new();
    let mut record_id = 0;

    let limit = if limit == 0 { usize::MAX } else { limit };

    while entries.len() < limit {
        let line_start = reader.virtual_position();
        line.clear();

        let bytes_read = reader.read_until(b'\n', &mut line)?;
        if bytes_read == 0 {
            break; // EOF
        }

        // Skip empty lines
        if line.len() <= 1 {
            continue;
        }

        if let Some((offset, len)) = tracepoint_offset_and_len(&line) {
            // Seek back to line start
            reader.seek(line_start)?;

            // Advance to tracepoint position
            if offset > 0 {
                io::copy(
                    &mut reader.by_ref().take(offset as u64),
                    &mut io::sink(),
                )?;
            }

            let tp_pos = reader.virtual_position();

            // Skip remaining bytes in line
            let remaining = line.len() as u64 - offset as u64;
            if remaining > 0 {
                io::copy(&mut reader.by_ref().take(remaining), &mut io::sink())?;
            }

            if len > 0 {
                entries.push(TracepointEntry {
                    record_id,
                    pos: tp_pos,
                    len,
                });
            }

            record_id += 1;
        }
    }

    Ok(entries)
}

/// Find the byte offset and length of the tracepoint tag (tp:Z:) in a PAF line
fn tracepoint_offset_and_len(line: &[u8]) -> Option<(usize, usize)> {
    const TAG: &[u8] = b"tp:Z:";
    let pos = line.windows(TAG.len()).position(|window| window == TAG)?;

    let start = pos + TAG.len();
    if start >= line.len() {
        return None;
    }

    // Count bytes until tab, newline, or carriage return
    let len = line[start..]
        .iter()
        .take_while(|&&b| b != b'\t' && b != b'\n' && b != b'\r')
        .count();

    if len == 0 {
        None
    } else {
        Some((start, len))
    }
}

/// Select random positions for testing (with fixed seed for reproducibility)
fn select_random_positions(entries: &[TracepointEntry], num_positions: usize) -> Vec<usize> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(42); // Fixed seed for reproducibility
    let mut indices: Vec<usize> = (0..entries.len()).collect();
    indices.shuffle(&mut rng);
    indices.truncate(num_positions.min(entries.len()));
    indices.sort_unstable(); // Sort for better cache behavior during reporting
    indices
}

fn benchmark_tracepoint_seeks(
    paf_path: &str,
    entries: &[TracepointEntry],
    test_positions: &[usize],
    iterations: usize,
) -> Result<SeekStats, Box<dyn std::error::Error>> {
    if test_positions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No positions available for benchmarking",
        )
        .into());
    }

    let file = File::open(paf_path)?;
    let mut reader = Reader::new(file);
    let mut buffer = Vec::new();
    let mut times = Vec::new();
    let mut successes = 0usize;
    let mut total_seeks = 0usize;

    for &idx in test_positions {
        let entry = entries[idx];
        if entry.len == 0 {
            continue;
        }

        buffer.resize(entry.len, 0);

        // Warm-up read
        if reader.seek(entry.pos).is_ok() && reader.read_exact(&mut buffer).is_ok() {
            // Benchmark iterations
            for _ in 0..iterations {
                let start = Instant::now();
                let seek_ok = reader.seek(entry.pos).is_ok();
                let read_ok = reader.read_exact(&mut buffer).is_ok();
                let elapsed = start.elapsed().as_micros() as f64;

                times.push(elapsed);
                total_seeks += 1;
                if seek_ok && read_ok {
                    successes += 1;
                }
            }
        }
    }

    if times.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No successful seeks during benchmarking",
        )
        .into());
    }

    // Calculate statistics
    let sum: f64 = times.iter().sum();
    let avg = sum / times.len() as f64;

    let variance: f64 = times.iter().map(|&t| (t - avg).powi(2)).sum::<f64>() / times.len() as f64;
    let stddev = variance.sqrt();

    // Sort for percentiles
    let mut sorted_times = times.clone();
    sorted_times.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let min = sorted_times[0];
    let max = sorted_times[sorted_times.len() - 1];
    let median = sorted_times[sorted_times.len() / 2];
    let p95 = sorted_times[(sorted_times.len() as f64 * 0.95) as usize];
    let p99 = sorted_times[(sorted_times.len() as f64 * 0.99) as usize];

    let success_rate = successes as f64 / total_seeks as f64;

    Ok(SeekStats {
        avg,
        stddev,
        median,
        min,
        max,
        p95,
        p99,
        success_rate,
    })
}

/// Save index to disk for faster subsequent runs
fn save_index(entries: &[TracepointEntry], path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Magic number and version
    writer.write_all(b"TPIDX")?;
    writer.write_all(&[1u8])?; // Version 1

    // Number of entries
    writer.write_all(&entries.len().to_le_bytes())?;

    // Write all entries
    for entry in entries {
        writer.write_all(&entry.to_bytes())?;
    }

    writer.flush()?;
    Ok(())
}

/// Load index from disk
fn load_index(path: &str) -> Result<Vec<TracepointEntry>, Box<dyn std::error::Error>> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 5];
    file.read_exact(&mut magic)?;

    if &magic != b"TPIDX" {
        return Err("Invalid index file: wrong magic number".into());
    }

    let mut version = [0u8; 1];
    file.read_exact(&mut version)?;
    if version[0] != 1 {
        return Err(format!("Unsupported index version: {}", version[0]).into());
    }

    let mut count_bytes = [0u8; 8];
    file.read_exact(&mut count_bytes)?;
    let count = usize::from_le_bytes(count_bytes);

    let mut entries = Vec::with_capacity(count);
    let mut entry_bytes = [0u8; 24];

    for _ in 0..count {
        file.read_exact(&mut entry_bytes)?;
        entries.push(TracepointEntry::from_bytes(&entry_bytes)?);
    }

    Ok(entries)
}
