use lib_bpaf::{
    is_binary_paf, BinaryPafHeader, StringTable, varint_size,
};
use std::env;
use std::fs::File;
use std::io::{self, BufReader, Seek};

/// Structure to hold detailed size analysis of a BPAF file
#[derive(Debug)]
struct BpafSizeAnalysis {
    // File sections in bytes (as stored in the file)
    total_file_size: u64,
    header_size: u64,
    string_table_size: u64,
    records_section_size: u64,

    // Metadata
    num_records: u64,
    num_strings: u64,
}

impl BpafSizeAnalysis {
    fn print_report(&self, header: &BinaryPafHeader) {
        // Calculate header field sizes
        let magic_size = 4u64;  // "BPAF" magic bytes
        let version_strategy_size = 2u64;  // version (1) + strategy (1)
        let num_records_varint = varint_size(self.num_records);
        let num_strings_varint = varint_size(self.num_strings);
        let tp_type_size = 1u64;  // TracepointType byte
        let complexity_metric_size = 1u64;  // ComplexityMetric byte
        let max_complexity_varint = varint_size(header.max_complexity());
        let distance_size = match header.distance() {
            lib_bpaf::Distance::Edit => 1u64,
            lib_bpaf::Distance::GapAffine { .. } => 1 + 3 * 4, // code + 3 i32 values
            lib_bpaf::Distance::GapAffine2p { .. } => 1 + 5 * 4, // code + 5 i32 values
        };

        println!("=== BPAF Header ==============");
        println!("  Magic bytes:                {:>4} bytes", magic_size);
        println!("  Version + strategy:         {:>4} bytes", version_strategy_size);
        println!("  Num records (varint):       {:>4} bytes  (value: {})", num_records_varint, self.num_records);
        println!("  Num strings (varint):       {:>4} bytes  (value: {})", num_strings_varint, self.num_strings);
        println!("  Tracepoint type:            {:>4} bytes  ({:?})", tp_type_size, header.tp_type());
        println!("  Complexity metric:          {:>4} bytes  ({:?})", complexity_metric_size, header.complexity_metric());
        println!("  Max complexity (varint):    {:>4} bytes  (value: {})", max_complexity_varint, header.max_complexity());
        println!("  Distance mode:              {:>4} bytes  ({:?})", distance_size, header.distance());
        println!("  Total:                      {:>4} bytes  ({:>6.2}%)",
            self.header_size,
            self.header_size as f64 / self.total_file_size as f64 * 100.0
        );

        println!("\n=== String Table ==============");
        println!("  Total size:                 {:>4} bytes  ({:>6.2}%)",
            self.string_table_size,
            self.string_table_size as f64 / self.total_file_size as f64 * 100.0
        );
        println!("  Avg bytes per string:       {:>7.2}",
            self.string_table_size as f64 / self.num_strings as f64
        );

        println!("\n=== Alignment Records ==============");
        println!("  Total size:                 {:>4} bytes  ({:>6.2}%)",
            self.records_section_size,
            self.records_section_size as f64 / self.total_file_size as f64 * 100.0
        );
        println!("  Avg bytes per record:       {:>7.2}",
            self.records_section_size as f64 / self.num_records as f64
        );

        println!("\n=== TOTAL ===");
        println!("Total file size:              {:>6} bytes  (100.00%)", self.total_file_size);
    }
}

fn analyze_bpaf_size(path: &str) -> io::Result<(BpafSizeAnalysis, BinaryPafHeader)> {
    let file = File::open(path)?;
    let total_file_size = file.metadata()?.len();
    let mut reader = BufReader::new(file);

    // Read header and track position
    let start_pos = reader.stream_position()?;
    let header = BinaryPafHeader::read(&mut reader)?;
    let header_end_pos = reader.stream_position()?;
    let header_size = header_end_pos - start_pos;

    // Read string table and track position
    let string_table_start = reader.stream_position()?;
    let _string_table = StringTable::read(&mut reader)?;
    let string_table_end = reader.stream_position()?;
    let string_table_size = string_table_end - string_table_start;

    // Calculate records section by subtraction
    let records_section_size = total_file_size - header_size - string_table_size;

    let analysis = BpafSizeAnalysis {
        total_file_size,
        header_size,
        string_table_size,
        records_section_size,
        num_records: header.num_records(),
        num_strings: header.num_strings(),
    };

    Ok((analysis, header))
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();

    if args.len() != 2 {
        eprintln!("Usage: {} <file.bpaf>", args[0]);
        eprintln!("\nAnalyzes a BPAF file and reports detailed size breakdown.");
        std::process::exit(1);
    }

    let path = &args[1];

    if !is_binary_paf(path)? {
        eprintln!("Error: '{}' is not a valid BPAF file", path);
        std::process::exit(1);
    }

    let (analysis, header) = analyze_bpaf_size(path)?;
    analysis.print_report(&header);

    Ok(())
}
