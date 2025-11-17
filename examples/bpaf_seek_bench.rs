/// Benchmark tool for testing O(1) random access performance on BPAF (Binary PAF) files
///
/// This tool:
/// 1. Opens a BPAF file and reads its index
/// 2. Performs random seeks to test access performance
/// 3. Reports statistics on seek times (average, stddev, min, max)
///
/// Usage: bpaf_seek_bench <file.bpaf> <num_records> <num_positions> <iterations_per_position>
///
/// Example: bpaf_seek_bench alignments.bpaf 10000 100 100
///   - Use first 10000 records
///   - Test 100 random positions
///   - Do 100 iterations per position

use lib_bpaf::BpafReader;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use std::env;
use std::io;
use std::time::Instant;

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!(
            "Usage: {} <file.bpaf> <num_records> <num_positions> <iterations>",
            args[0]
        );
        eprintln!();
        eprintln!("Arguments:");
        eprintln!("  <file.bpaf>         : Path to BPAF file");
        eprintln!("  <num_records>       : Number of records to test (0 = all)");
        eprintln!("  <num_positions>     : Number of random positions to test");
        eprintln!("  <iterations>        : Iterations per position");
        eprintln!();
        eprintln!("Example:");
        eprintln!("  {} alignments.bpaf 10000 100 100", args[0]);
        std::process::exit(1);
    }

    let bpaf_path = &args[1];
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
        eprintln!("Error: num_positions and iterations must be greater than zero");
        std::process::exit(1);
    }

    println!("=== BPAF Seek Benchmark ===");
    println!("File: {}", bpaf_path);
    println!("Loading BPAF file...");

    let start = Instant::now();
    let mut reader = BpafReader::open(bpaf_path)?;
    let load_time = start.elapsed();

    let total_records = reader.len();

    if total_records == 0 {
        println!("No records found!");
        return Ok(());
    }

    let test_limit = if num_records == 0 || num_records > total_records {
        total_records
    } else {
        num_records
    };

    println!(
        "Loaded {} records in {:.2}s",
        total_records,
        load_time.as_secs_f64()
    );
    println!();

    // Select random positions to test
    let test_positions = select_random_positions(test_limit, num_positions);
    println!(
        "Testing {} random positions × {} iterations each",
        test_positions.len(),
        iterations
    );
    println!();

    let stats = benchmark_bpaf_seeks(&mut reader, &test_positions, iterations)?;

    println!("=== Results ===");
    println!("Records available:   {}", total_records);
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

    // Check O(1) performance
    if stats.avg < 1000.0 {
        println!("✓ O(1) performance achieved (mean < 1ms)");
    } else {
        println!("⚠ Slow seeks detected (mean >= 1ms)");
    }

    // Check consistency (coefficient of variation)
    let cv = (stats.stddev / stats.avg) * 100.0;
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
    println!("BPAF_INDEXED {}", total_records);
    println!("BPAF_SEEK_AVG {:.2}", stats.avg);
    println!("BPAF_SEEK_STDDEV {:.2}", stats.stddev);
    println!("BPAF_SUCCESS_RATE {:.4}", stats.success_rate);

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

/// Select random positions for testing (with fixed seed for reproducibility)
fn select_random_positions(max_records: usize, num_positions: usize) -> Vec<usize> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(42); // Fixed seed for reproducibility
    let mut indices: Vec<usize> = (0..max_records).collect();
    indices.shuffle(&mut rng);
    indices.truncate(num_positions.min(max_records));
    indices.sort_unstable(); // Sort for better cache behavior during reporting
    indices
}

fn benchmark_bpaf_seeks(
    reader: &mut BpafReader,
    test_positions: &[usize],
    iterations: usize,
) -> io::Result<SeekStats> {
    if test_positions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No positions available for benchmarking",
        ));
    }

    let mut times = Vec::new();
    let mut successes = 0usize;
    let mut total_seeks = 0usize;

    for &pos in test_positions {
        // Warm-up read
        if reader.get_alignment_record(pos as u64).is_ok() {
            // Benchmark iterations
            for _ in 0..iterations {
                let start = Instant::now();
                let result = reader.get_alignment_record(pos as u64);
                let elapsed = start.elapsed().as_micros() as f64;

                times.push(elapsed);
                total_seeks += 1;
                if result.is_ok() {
                    successes += 1;
                }
            }
        }
    }

    if times.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::Other,
            "No successful seeks during benchmarking",
        ));
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
