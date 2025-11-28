use lib_wfa2::affine_wavefront::Distance;
use tpa::{compress_paf_to_tpa, CompressionConfig, CompressionStrategy};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::init();
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: auto_check <input.paf> <output.tpa>");
        std::process::exit(1);
    }

    compress_paf_to_tpa(
        &args[1],
        &args[2],
        CompressionConfig::new()
            .strategy(CompressionStrategy::Automatic(3, 0)) // 0 = analyze entire file
            .distance(Distance::GapAffine {
                mismatch: 8,
                gap_opening: 5,
                gap_extension: 2,
            }),
    )?;

    Ok(())
}
