use std::env;
use lib_bpaf::{compress_paf_with_tracepoints_dual, CompressionStrategy, TracepointType, ComplexityMetric, Distance};

fn parse_strategy(s: &str) -> Result<CompressionStrategy, String> {
    let parts: Vec<&str> = s.split(',').collect();
    let name = parts[0];
    let level = if parts.len() > 1 {
        parts[1].parse::<i32>().map_err(|_| format!("Invalid zstd level: {}", parts[1]))?
    } else {
        3 // default
    };

    match name {
        "raw" => Ok(CompressionStrategy::Raw(level)),
        "zigzag-delta" => Ok(CompressionStrategy::ZigzagDelta(level)),
        "2d-delta" => Ok(CompressionStrategy::TwoDimDelta(level)),
        "rle" => Ok(CompressionStrategy::RunLength(level)),
        "bit-packed" => Ok(CompressionStrategy::BitPacked(level)),
        "delta-of-delta" => Ok(CompressionStrategy::DeltaOfDelta(level)),
        "frame-of-reference" => Ok(CompressionStrategy::FrameOfReference(level)),
        "hybrid-rle" => Ok(CompressionStrategy::HybridRLE(level)),
        "offset-joint" => Ok(CompressionStrategy::OffsetJoint(level)),
        "xor-delta" => Ok(CompressionStrategy::XORDelta(level)),
        "dictionary" => Ok(CompressionStrategy::Dictionary(level)),
        "simple8" => Ok(CompressionStrategy::Simple8(level)),
        "stream-vbyte" => Ok(CompressionStrategy::StreamVByte(level)),
        "fastpfor" => Ok(CompressionStrategy::FastPFOR(level)),
        "cascaded" => Ok(CompressionStrategy::Cascaded(level)),
        "simple8b-full" => Ok(CompressionStrategy::Simple8bFull(level)),
        "selective-rle" => Ok(CompressionStrategy::SelectiveRLE(level)),
        _ => Err(format!("Unknown strategy: {}", name)),
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 4 {
        eprintln!("Usage: {} <input.paf> <output.bpaf> <first_strategy> <second_strategy> [tp_type] [max_complexity] [metric]", args[0]);
        eprintln!();
        eprintln!("Examples:");
        eprintln!("  {} input.paf output.bpaf raw zigzag-delta", args[0]);
        eprintln!("  {} input.paf output.bpaf raw,3 zigzag-delta,6", args[0]);
        eprintln!("  {} input.paf output.bpaf 2d-delta 2d-delta standard 32 edit-distance", args[0]);
        eprintln!();
        eprintln!("Strategies: raw, zigzag-delta, 2d-delta, rle, bit-packed, delta-of-delta,");
        eprintln!("           frame-of-reference, hybrid-rle, offset-joint, xor-delta, dictionary,");
        eprintln!("           simple8, stream-vbyte, fastpfor, cascaded, simple8b-full, selective-rle");
        std::process::exit(1);
    }

    let input_path = &args[1];
    let output_path = &args[2];
    let first_strategy = parse_strategy(&args[3])?;
    let second_strategy = parse_strategy(&args[4])?;

    let tp_type = if args.len() > 5 {
        match args[5].as_str() {
            "standard" => TracepointType::Standard,
            "variable" => TracepointType::Variable,
            "mixed" => TracepointType::Mixed,
            "fastga" => TracepointType::Fastga,
            _ => return Err(format!("Unknown tracepoint type: {}", args[5]).into()),
        }
    } else {
        TracepointType::Standard
    };

    let max_complexity = if args.len() > 6 {
        args[6].parse::<u64>()?
    } else {
        32
    };

    let complexity_metric = if args.len() > 7 {
        match args[7].as_str() {
            "edit-distance" => ComplexityMetric::EditDistance,
            "diagonal-distance" => ComplexityMetric::DiagonalDistance,
            _ => return Err(format!("Unknown metric: {}", args[7]).into()),
        }
    } else {
        ComplexityMetric::EditDistance
    };

    println!("Compressing {} to {}", input_path, output_path);
    println!("First strategy:  {:?}", first_strategy);
    println!("Second strategy: {:?}", second_strategy);
    println!("Tracepoint type: {:?}", tp_type);
    println!("Max complexity:  {}", max_complexity);
    println!("Metric:          {:?}", complexity_metric);

    compress_paf_with_tracepoints_dual(
        input_path,
        output_path,
        first_strategy,
        second_strategy,
        tp_type,
        max_complexity,
        complexity_metric,
        Distance::Edit,
    )?;

    println!("✓ Compression complete");
    Ok(())
}
