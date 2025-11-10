use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use lib_bpaf::BpafReader;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <original_paf> <bpaf_file> <record_id>", args[0]);
        std::process::exit(1);
    }

    let paf_path = &args[1];
    let bpaf_path = &args[2];
    let record_id: u64 = args[3].parse().expect("Invalid record_id");

    // Read original PAF to get expected tracepoints
    let file = File::open(paf_path).expect("Failed to open PAF");
    let reader = BufReader::new(file);

    let mut current_line = 0u64;
    let mut expected_tp: Option<String> = None;

    for line in reader.lines() {
        if current_line == record_id {
            let line = line.expect("Failed to read line");
            // Extract tp:Z: tag
            for field in line.split('\t') {
                if field.starts_with("tp:Z:") {
                    expected_tp = Some(field.strip_prefix("tp:Z:").unwrap().to_string());
                    break;
                }
            }
            break;
        }
        current_line += 1;
    }

    let expected = expected_tp.expect("Record not found in PAF or no tp:Z: tag");

    // Read from BPAF
    let mut bpaf_reader = BpafReader::open(bpaf_path).expect("Failed to open BPAF");
    let (tracepoints, _tp_type, _complexity, _max) = bpaf_reader
        .get_tracepoints(record_id)
        .expect("Failed to get tracepoints from BPAF");

    // Convert tracepoints to string format
    let actual = match tracepoints {
        lib_bpaf::TracepointData::Standard(tps) | lib_bpaf::TracepointData::Fastga(tps) => {
            tps.iter()
                .map(|(a, b)| format!("{},{}", a, b))
                .collect::<Vec<_>>()
                .join(";")
        }
        lib_bpaf::TracepointData::Variable(tps) => {
            tps.iter()
                .map(|(a, b)| {
                    if let Some(b_val) = b {
                        format!("{},{}", a, b_val)
                    } else {
                        format!("{}", a)
                    }
                })
                .collect::<Vec<_>>()
                .join(";")
        }
        lib_bpaf::TracepointData::Mixed(_) => {
            panic!("Mixed tracepoints not expected in test data");
        }
    };

    // Compare
    if expected == actual {
        println!("✓ Record {} tracepoints verified", record_id);
        println!("  Values: {}", if actual.len() > 80 {
            format!("{}... ({} chars)", &actual[..80], actual.len())
        } else {
            actual
        });
    } else {
        println!("✗ Record {} tracepoints MISMATCH!", record_id);
        println!("  Expected: {}", if expected.len() > 80 {
            format!("{}...", &expected[..80])
        } else {
            expected.clone()
        });
        println!("  Actual:   {}", if actual.len() > 80 {
            format!("{}...", &actual[..80])
        } else {
            actual.clone()
        });
        std::process::exit(1);
    }
}
