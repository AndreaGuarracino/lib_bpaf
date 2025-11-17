use std::env;
use std::fs::File;
use lib_bpaf::BpafReader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: {} <bpaf_file>", args[0]);
        std::process::exit(1);
    }

    let bpaf_path = &args[1];
    let mut reader = BpafReader::open(bpaf_path)?;

    let header = reader.header();

    // Print both strategies (extract from Ok() wrapper and remove Debug formatting)
    let first_strat = header.first_strategy()?;
    let second_strat = header.second_strategy()?;

    // Convert to simple string representation (remove the compression level number)
    let first_name = format!("{:?}", first_strat).split('(').next().unwrap_or("Unknown").to_lowercase();
    let second_name = format!("{:?}", second_strat).split('(').next().unwrap_or("Unknown").to_lowercase();

    println!("{}\t{}", first_name, second_name);

    Ok(())
}
