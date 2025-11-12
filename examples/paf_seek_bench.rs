use noodles::bgzf::io::Reader;
use noodles::bgzf::VirtualPosition;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, Read};
use std::time::Instant;

#[derive(Clone, Copy)]
struct CigarEntry {
    pos: VirtualPosition,
    len: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        eprintln!(
            "Usage: {} <paf.bgz> <num_records> <stride> <iterations>",
            args[0]
        );
        std::process::exit(1);
    }

    let paf_path = &args[1];
    let num_records: usize = args[2]
        .parse()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid num_records"))?;
    let stride: usize = args[3]
        .parse()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid stride"))?;
    let iterations: usize = args[4]
        .parse()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid iterations"))?;

    if stride == 0 || iterations == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Stride and iterations must be greater than zero",
        )
        .into());
    }

    let entries = build_cigar_index(paf_path, num_records)?;

    if entries.is_empty() {
        println!("PAF_INDEXED 0");
        println!("PAF_SEEK N/A");
        return Ok(());
    }

    let avg_seek = benchmark_cigar_seeks(paf_path, &entries, stride, iterations)?;

    println!("PAF_INDEXED {}", entries.len());
    println!("PAF_SEEK {:.2}", avg_seek);

    Ok(())
}

fn build_cigar_index(
    paf_path: &str,
    limit: usize,
) -> Result<Vec<CigarEntry>, Box<dyn std::error::Error>> {
    let file = File::open(paf_path)?;
    let mut reader = Reader::new(file);
    let mut line = Vec::new();
    let mut entries = Vec::new();

    while entries.len() < limit {
        let line_start = reader.virtual_position();
        line.clear();

        let bytes_read = reader.read_until(b'\n', &mut line)?;
        if bytes_read == 0 {
            break;
        }

        if let Some((offset, len)) = cigar_offset_and_len(&line) {
            reader.seek(line_start)?;
            if offset > 0 {
                io::copy(
                    &mut reader.by_ref().take(offset as u64),
                    &mut io::sink(),
                )?;
            }
            let cigar_pos = reader.virtual_position();

            let remaining = line.len() as u64 - offset as u64;
            if remaining > 0 {
                io::copy(&mut reader.by_ref().take(remaining), &mut io::sink())?;
            }

            if len > 0 {
                entries.push(CigarEntry {
                    pos: cigar_pos,
                    len,
                });
            }
        }
    }

    Ok(entries)
}

fn cigar_offset_and_len(line: &[u8]) -> Option<(usize, usize)> {
    const TAG: &[u8] = b"cg:Z:";
    let pos = line.windows(TAG.len()).position(|window| window == TAG)?;

    let start = pos + TAG.len();
    if start >= line.len() {
        return None;
    }

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

fn benchmark_cigar_seeks(
    paf_path: &str,
    entries: &[CigarEntry],
    stride: usize,
    iterations: usize,
) -> Result<f64, Box<dyn std::error::Error>> {
    if entries.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No entries available for benchmarking",
        )
        .into());
    }

    let file = File::open(paf_path)?;
    let mut reader = Reader::new(file);
    let mut buffer = Vec::new();
    let mut total = 0.0;
    let mut samples = 0usize;

    for idx in (0..entries.len()).step_by(stride) {
        let entry = entries[idx];
        if entry.len == 0 {
            continue;
        }

        buffer.resize(entry.len, 0);

        // Warm-up read
        reader.seek(entry.pos)?;
        reader.read_exact(&mut buffer)?;

        let mut subtotal = 0.0;
        for _ in 0..iterations {
            let start = Instant::now();
            reader.seek(entry.pos)?;
            reader.read_exact(&mut buffer)?;
            subtotal += start.elapsed().as_micros() as f64;
        }

        total += subtotal / iterations as f64;
        samples += 1;
    }

    if samples == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "No samples collected during benchmarking",
        )
        .into());
    }

    Ok(total / samples as f64)
}
