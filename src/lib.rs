/// Binary PAF format for efficient storage of sequence alignments with tracepoints
///
/// Format: [Header] → [Records] → [StringTable]
/// - Header: Magic "BPAF" + version + metadata
/// - Records: Core PAF fields + compressed tracepoints
/// - StringTable: Deduplicated sequence names with lengths
///
/// Compression strategy:
/// - Delta encoding: positions stored as deltas from previous position
/// - Varint encoding: variable-length integer encoding
/// - Zstd level 3: fast compression with good ratios
use flate2::read::MultiGzDecoder;
use lib_tracepoints::{
    cigar_to_mixed_tracepoints, cigar_to_tracepoints, cigar_to_tracepoints_fastga,
    cigar_to_variable_tracepoints, ComplexityMetric, MixedRepresentation, TracepointType,
};
use log::{debug, error, info};
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

// ============================================================================
// COMPRESSION STRATEGY
// ============================================================================

/// Compression strategy for binary PAF format
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressionStrategy {
    /// Delta (FastGA-aware) + Varint + Zstd level 3
    /// - FastGA: raw num_differences values (naturally small)
    /// - Standard: delta-encoded query_bases
    Varint,
    /// No delta encoding + Varint + Zstd level 3
    /// - All types: raw first values (test if delta hurts compression)
    VarintRaw,
    /// Delta (FastGA-aware) + Huffman + Zstd
    /// - FastGA: raw num_differences (exploits natural skew)
    /// - Standard: delta-encoded query_bases
    Huffman,
    /// No delta encoding + Huffman + Zstd
    /// - All types: raw first values (maximum entropy coding)
    HuffmanRaw,
    /// Hybrid: Delta + Varint + Byte-Huffman + Zstd
    /// - Encodes integers with varint (produces bytes)
    /// - Then applies Huffman on byte distribution
    /// - Tests if byte-level Huffman helps Zstd
    VarintHuffman,
    /// Hybrid: Raw + Varint + Byte-Huffman + Zstd
    /// - No delta, varint on raw values
    /// - Then byte-level Huffman
    VarintHuffmanRaw,
    /// Analyze data distribution and auto-select best strategy
    /// - Measures top 3 symbol coverage
    /// - High skew (>60%): Huffman
    /// - Flat distribution: Varint
    /// - Compares delta vs raw for Huffman decision
    Smart,
}

impl CompressionStrategy {
    /// Parse strategy from string
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s.to_lowercase().as_str() {
            "varint" => Ok(CompressionStrategy::Varint),
            "varint-raw" => Ok(CompressionStrategy::VarintRaw),
            "huffman" => Ok(CompressionStrategy::Huffman),
            "huffman-raw" => Ok(CompressionStrategy::HuffmanRaw),
            "varint-huffman" => Ok(CompressionStrategy::VarintHuffman),
            "varint-huffman-raw" => Ok(CompressionStrategy::VarintHuffmanRaw),
            "smart" => Ok(CompressionStrategy::Smart),
            _ => Err(format!("Unknown compression strategy: {}", s)),
        }
    }

    /// Get all available strategies
    pub fn variants() -> &'static [&'static str] {
        &["varint", "varint-raw", "huffman", "huffman-raw", "varint-huffman", "varint-huffman-raw", "smart"]
    }

    /// Convert to strategy code for file header
    fn to_code(&self) -> u8 {
        match self {
            CompressionStrategy::Varint => 0,
            CompressionStrategy::VarintRaw => 1,
            CompressionStrategy::Huffman => 2,
            CompressionStrategy::HuffmanRaw => 3,
            CompressionStrategy::VarintHuffman => 4,
            CompressionStrategy::VarintHuffmanRaw => 5,
            CompressionStrategy::Smart => unreachable!("Smart strategy should be resolved before writing"),
        }
    }

    /// Parse from strategy code
    fn from_code(code: u8) -> io::Result<Self> {
        match code {
            0 => Ok(CompressionStrategy::Varint),
            1 => Ok(CompressionStrategy::VarintRaw),
            2 => Ok(CompressionStrategy::Huffman),
            3 => Ok(CompressionStrategy::HuffmanRaw),
            4 => Ok(CompressionStrategy::VarintHuffman),
            5 => Ok(CompressionStrategy::VarintHuffmanRaw),
            _ => Err(io::Error::new(io::ErrorKind::InvalidData, format!("Unknown strategy code: {}", code))),
        }
    }

    /// Check if strategy requires Huffman codecs
    fn requires_codecs(&self) -> bool {
        matches!(self,
            CompressionStrategy::Huffman |
            CompressionStrategy::HuffmanRaw |
            CompressionStrategy::VarintHuffman |
            CompressionStrategy::VarintHuffmanRaw
        )
    }

    /// Check if strategy uses delta encoding
    fn uses_delta(&self) -> bool {
        matches!(self,
            CompressionStrategy::Varint |
            CompressionStrategy::Huffman |
            CompressionStrategy::VarintHuffman
        )
    }
}

impl std::fmt::Display for CompressionStrategy {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CompressionStrategy::Varint => write!(f, "varint"),
            CompressionStrategy::VarintRaw => write!(f, "varint-raw"),
            CompressionStrategy::Huffman => write!(f, "huffman"),
            CompressionStrategy::HuffmanRaw => write!(f, "huffman-raw"),
            CompressionStrategy::VarintHuffman => write!(f, "varint-huffman"),
            CompressionStrategy::VarintHuffmanRaw => write!(f, "varint-huffman-raw"),
            CompressionStrategy::Smart => write!(f, "smart"),
        }
    }
}

// ============================================================================
// VARINT ENCODING (LEB128)
// ============================================================================

/// Encode an unsigned integer as a varint
fn encode_varint(mut value: u64) -> Vec<u8> {
    let mut bytes = Vec::new();
    loop {
        let mut byte = (value & 0x7F) as u8;
        value >>= 7;
        if value != 0 {
            byte |= 0x80; // Set continuation bit
        }
        bytes.push(byte);
        if value == 0 {
            break;
        }
    }
    bytes
}

/// Write a varint to a writer
fn write_varint<W: Write>(writer: &mut W, value: u64) -> io::Result<usize> {
    let bytes = encode_varint(value);
    writer.write_all(&bytes)?;
    Ok(bytes.len())
}

/// Read a varint from a reader
fn read_varint<R: Read>(reader: &mut R) -> io::Result<u64> {
    let mut value: u64 = 0;
    let mut shift = 0;
    loop {
        let mut byte_buf = [0u8; 1];
        reader.read_exact(&mut byte_buf)?;
        let byte = byte_buf[0];
        value |= ((byte & 0x7F) as u64) << shift;
        if byte & 0x80 == 0 {
            break;
        }
        shift += 7;
        if shift >= 64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Varint too long",
            ));
        }
    }
    Ok(value)
}

// ============================================================================
// HUFFMAN CODING
// ============================================================================

/// Huffman tree node
#[derive(Clone)]
enum HuffmanNode {
    Leaf { value: i64, frequency: u64 },
    Internal { left: Box<HuffmanNode>, right: Box<HuffmanNode>, frequency: u64 },
}

impl HuffmanNode {
    fn frequency(&self) -> u64 {
        match self {
            Self::Leaf { frequency, .. } => *frequency,
            Self::Internal { frequency, .. } => *frequency,
        }
    }
}

impl Eq for HuffmanNode {}
impl PartialEq for HuffmanNode {
    fn eq(&self, other: &Self) -> bool {
        self.frequency() == other.frequency()
    }
}
impl Ord for HuffmanNode {
    fn cmp(&self, other: &Self) -> Ordering {
        other.frequency().cmp(&self.frequency()) // Reverse for min-heap
    }
}
impl PartialOrd for HuffmanNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Variable-length bit string (for backwards compatibility reading)
#[derive(Clone, Debug)]
struct BitString {
    bits: Vec<bool>,
}

/// Huffman tree for decoding (backwards compatibility only)
#[derive(Clone)]
struct HuffmanTree {
    root: HuffmanNode,
    encode_table: HashMap<i64, BitString>,
}

impl HuffmanTree {
    fn decode(&self, bits: &mut BitReader) -> io::Result<i64> {
        let mut node = &self.root;
        loop {
            match node {
                HuffmanNode::Leaf { value, .. } => return Ok(*value),
                HuffmanNode::Internal { left, right, .. } => {
                    let bit = bits.read_bit()?;
                    node = if bit { right } else { left };
                }
            }
        }
    }

    fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let num_symbols = read_varint(reader)? as usize;
        let mut encode_table = HashMap::new();

        for _ in 0..num_symbols {
            let zigzag = read_varint(reader)?;
            let value = ((zigzag >> 1) as i64) ^ -((zigzag & 1) as i64);

            let mut len_buf = [0u8; 1];
            reader.read_exact(&mut len_buf)?;
            let bit_len = len_buf[0] as usize;

            let num_bytes = (bit_len + 7) / 8;
            let mut bytes = vec![0u8; num_bytes];
            reader.read_exact(&mut bytes)?;

            let mut bits = Vec::with_capacity(bit_len);
            for i in 0..bit_len {
                let byte_idx = i / 8;
                let bit_idx = 7 - (i % 8);
                bits.push((bytes[byte_idx] >> bit_idx) & 1 == 1);
            }

            encode_table.insert(value, BitString { bits });
        }

        let root = rebuild_tree(&encode_table);
        Ok(Self { root, encode_table })
    }

    fn num_symbols(&self) -> usize {
        self.encode_table.len()
    }
}

fn rebuild_tree(encode_table: &HashMap<i64, BitString>) -> HuffmanNode {
    if encode_table.len() == 1 {
        let (&value, _) = encode_table.iter().next().unwrap();
        return HuffmanNode::Leaf { value, frequency: 1 };
    }

    let mut root = HuffmanNode::Internal {
        frequency: 0,
        left: Box::new(HuffmanNode::Leaf { value: 0, frequency: 0 }),
        right: Box::new(HuffmanNode::Leaf { value: 0, frequency: 0 }),
    };

    for (&value, code) in encode_table {
        let mut node = &mut root;
        for (i, &bit) in code.bits.iter().enumerate() {
            let is_last = i == code.bits.len() - 1;
            match node {
                HuffmanNode::Internal { left, right, .. } => {
                    let child = if bit { right.as_mut() } else { left.as_mut() };
                    if is_last {
                        *child = HuffmanNode::Leaf { value, frequency: 1 };
                        break;
                    } else {
                        if matches!(child, HuffmanNode::Leaf { .. }) {
                            *child = HuffmanNode::Internal {
                                frequency: 0,
                                left: Box::new(HuffmanNode::Leaf { value: 0, frequency: 0 }),
                                right: Box::new(HuffmanNode::Leaf { value: 0, frequency: 0 }),
                            };
                        }
                        node = child;
                    }
                }
                _ => panic!("Unexpected leaf node in tree traversal"),
            }
        }
    }
    root
}

struct BitReader<'a> {
    bytes: &'a [u8],
    byte_pos: usize,
    bit_pos: u8,
}

impl<'a> BitReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, byte_pos: 0, bit_pos: 0 }
    }

    fn read_bit(&mut self) -> io::Result<bool> {
        if self.byte_pos >= self.bytes.len() {
            return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "BitReader: end of stream"));
        }
        let bit = (self.bytes[self.byte_pos] >> (7 - self.bit_pos)) & 1 == 1;
        self.bit_pos += 1;
        if self.bit_pos == 8 {
            self.byte_pos += 1;
            self.bit_pos = 0;
        }
        Ok(bit)
    }
}

impl HuffmanTree {
    /// Build Huffman tree from value frequencies
    fn build(frequencies: &HashMap<i64, u64>) -> Result<Self, String> {
        if frequencies.is_empty() {
            return Err("Cannot build Huffman tree from empty frequency map".to_string());
        }

        use std::collections::BinaryHeap;

        let mut heap: BinaryHeap<HuffmanNode> = frequencies
            .iter()
            .map(|(&value, &frequency)| HuffmanNode::Leaf { value, frequency })
            .collect();

        while heap.len() > 1 {
            let left = heap.pop().unwrap();
            let right = heap.pop().unwrap();
            let combined_freq = left.frequency() + right.frequency();
            heap.push(HuffmanNode::Internal {
                left: Box::new(left),
                right: Box::new(right),
                frequency: combined_freq,
            });
        }

        let root = heap.pop().unwrap();
        let encode_table = Self::generate_code_table(&root);

        // Debug: log code length statistics
        let code_lengths: Vec<usize> = encode_table.values().map(|bs| bs.bits.len()).collect();
        if !code_lengths.is_empty() {
            let min_len = *code_lengths.iter().min().unwrap();
            let max_len = *code_lengths.iter().max().unwrap();
            let avg_len = code_lengths.iter().sum::<usize>() as f64 / code_lengths.len() as f64;
            debug!("  Code lengths: min={}, max={}, avg={:.2}", min_len, max_len, avg_len);
        }

        Ok(Self { root, encode_table })
    }

    /// Generate code table from Huffman tree
    fn generate_code_table(node: &HuffmanNode) -> HashMap<i64, BitString> {
        let mut table = HashMap::new();
        Self::generate_codes_recursive(node, Vec::new(), &mut table);
        table
    }

    fn generate_codes_recursive(
        node: &HuffmanNode,
        current_code: Vec<bool>,
        table: &mut HashMap<i64, BitString>,
    ) {
        match node {
            HuffmanNode::Leaf { value, .. } => {
                let bits = if current_code.is_empty() {
                    vec![false] // Single symbol gets code "0"
                } else {
                    current_code
                };
                table.insert(*value, BitString { bits });
            }
            HuffmanNode::Internal { left, right, .. } => {
                let mut left_code = current_code.clone();
                left_code.push(false);
                Self::generate_codes_recursive(left, left_code, table);

                let mut right_code = current_code;
                right_code.push(true);
                Self::generate_codes_recursive(right, right_code, table);
            }
        }
    }

    /// Encode a value to bit sequence
    fn encode(&self, value: i64) -> Result<&BitString, String> {
        self.encode_table
            .get(&value)
            .ok_or_else(|| format!("Value {} not in Huffman code table", value))
    }

    /// Write Huffman tree to file
    fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        write_varint(writer, self.encode_table.len() as u64)?;

        for (&value, code) in &self.encode_table {
            // Zigzag encode signed value
            let zigzag = ((value << 1) ^ (value >> 63)) as u64;
            write_varint(writer, zigzag)?;

            // Write bit length
            writer.write_all(&[code.bits.len() as u8])?;

            // Pack bits into bytes
            let num_bytes = (code.bits.len() + 7) / 8;
            let mut bytes = vec![0u8; num_bytes];
            for (i, &bit) in code.bits.iter().enumerate() {
                if bit {
                    bytes[i / 8] |= 1 << (7 - (i % 8));
                }
            }
            writer.write_all(&bytes)?;
        }

        Ok(())
    }
}

struct BitWriter {
    bytes: Vec<u8>,
    current_byte: u8,
    bit_pos: u8,
}

impl BitWriter {
    fn new() -> Self {
        Self {
            bytes: Vec::new(),
            current_byte: 0,
            bit_pos: 0,
        }
    }

    fn write_bit(&mut self, bit: bool) {
        if bit {
            self.current_byte |= 1 << (7 - self.bit_pos);
        }
        self.bit_pos += 1;
        if self.bit_pos == 8 {
            self.bytes.push(self.current_byte);
            self.current_byte = 0;
            self.bit_pos = 0;
        }
    }

    fn write_bits(&mut self, bits: &[bool]) {
        for &bit in bits {
            self.write_bit(bit);
        }
    }

    fn finish(mut self) -> Vec<u8> {
        if self.bit_pos > 0 {
            self.bytes.push(self.current_byte);
        }
        self.bytes
    }
}

/// Adaptive codec for encoding/decoding with Huffman
#[derive(Clone)]
struct AdaptiveCodec {
    huffman_tree: Option<HuffmanTree>,
}

impl AdaptiveCodec {
    /// Build codec from value frequencies
    fn build(frequencies: &HashMap<i64, u64>) -> Result<Self, String> {
        if frequencies.is_empty() {
            return Err("Cannot build codec from empty frequency map".to_string());
        }
        let huffman_tree = HuffmanTree::build(frequencies)?;
        Ok(Self {
            huffman_tree: Some(huffman_tree),
        })
    }

    /// Encode a value
    fn encode(&self, value: i64) -> Result<&BitString, String> {
        self.huffman_tree.as_ref().unwrap().encode(value)
    }

    fn decode(&self, bits: &mut BitReader) -> io::Result<i64> {
        self.huffman_tree.as_ref().unwrap().decode(bits)
    }

    fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let huffman_tree = HuffmanTree::read(reader)?;
        Ok(Self {
            huffman_tree: Some(huffman_tree),
        })
    }

    fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        self.huffman_tree.as_ref().unwrap().write(writer)
    }

    fn num_symbols(&self) -> usize {
        self.huffman_tree
            .as_ref()
            .map(|t| t.num_symbols())
            .unwrap_or(0)
    }
}

// ============================================================================
// BINARY PAF FORMAT
// ============================================================================

const BINARY_MAGIC: &[u8; 4] = b"BPAF";
const FLAG_COMPRESSED: u8 = 0x01;
const FLAG_ADAPTIVE: u8 = 0x02;

/// Binary PAF header
#[derive(Debug)]
pub struct BinaryPafHeader {
    version: u8,
    flags: u8,
    num_records: u64,
    num_strings: u64,
}

impl BinaryPafHeader {
    fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(BINARY_MAGIC)?;
        writer.write_all(&[self.version, self.flags])?;
        write_varint(writer, self.num_records)?;
        write_varint(writer, self.num_strings)?;
        Ok(())
    }

    fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != BINARY_MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid magic"));
        }
        let mut ver_flags = [0u8; 2];
        reader.read_exact(&mut ver_flags)?;
        Ok(Self {
            version: ver_flags[0],
            flags: ver_flags[1],
            num_records: read_varint(reader)?,
            num_strings: read_varint(reader)?,
        })
    }
}

/// String table for deduplicating sequence names
#[derive(Debug, Default)]
pub struct StringTable {
    strings: Vec<String>,
    lengths: Vec<u64>,
    index: HashMap<String, u64>,
}

impl StringTable {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn intern(&mut self, s: &str, length: u64) -> u64 {
        if let Some(&id) = self.index.get(s) {
            id
        } else {
            let id = self.strings.len() as u64;
            self.strings.push(s.to_string());
            self.lengths.push(length);
            self.index.insert(s.to_string(), id);
            id
        }
    }

    pub fn get(&self, id: u64) -> Option<&str> {
        self.strings.get(id as usize).map(|s| s.as_str())
    }

    pub fn get_length(&self, id: u64) -> Option<u64> {
        self.lengths.get(id as usize).copied()
    }

    pub fn len(&self) -> usize {
        self.strings.len()
    }

    pub fn is_empty(&self) -> bool {
        self.strings.is_empty()
    }

    fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        write_varint(writer, self.strings.len() as u64)?;
        for (s, &len) in self.strings.iter().zip(self.lengths.iter()) {
            write_varint(writer, s.len() as u64)?;
            writer.write_all(s.as_bytes())?;
            write_varint(writer, len)?;
        }
        Ok(())
    }

    fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let num_strings = read_varint(reader)? as usize;
        let mut strings = Vec::with_capacity(num_strings);
        let mut lengths = Vec::with_capacity(num_strings);
        let mut index = HashMap::new();

        for id in 0..num_strings {
            let len = read_varint(reader)? as usize;
            let mut buf = vec![0u8; len];
            reader.read_exact(&mut buf)?;
            let s = String::from_utf8(buf)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let seq_len = read_varint(reader)?;
            index.insert(s.clone(), id as u64);
            strings.push(s);
            lengths.push(seq_len);
        }
        Ok(Self {
            strings,
            lengths,
            index,
        })
    }
}

pub struct AlignmentRecord {
    pub query_name_id: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub strand: char,
    pub target_name_id: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub residue_matches: u64,
    pub alignment_block_len: u64,
    pub mapping_quality: u8,
    pub tp_type: TracepointType,
    pub complexity_metric: ComplexityMetric,
    pub max_complexity: u64,
    pub tracepoints: TracepointData,
    pub tags: Vec<Tag>,
}

#[derive(Debug, Clone)]
pub enum TracepointData {
    Standard(Vec<(u64, u64)>),
    Mixed(Vec<MixedTracepointItem>),
    Variable(Vec<(u64, Option<u64>)>),
    Fastga(Vec<(u64, u64)>),
}

#[derive(Debug, Clone)]
pub enum MixedTracepointItem {
    Tracepoint(u64, u64),
    CigarOp(u64, char),
}

impl From<&MixedRepresentation> for MixedTracepointItem {
    fn from(mr: &MixedRepresentation) -> Self {
        match mr {
            MixedRepresentation::Tracepoint(a, b) => Self::Tracepoint(*a as u64, *b as u64),
            MixedRepresentation::CigarOp(len, op) => Self::CigarOp(*len as u64, *op),
        }
    }
}

#[derive(Debug)]
pub struct Tag {
    pub key: [u8; 2],
    pub tag_type: u8,
    pub value: TagValue,
}

#[derive(Debug)]
pub enum TagValue {
    Int(i64),
    Float(f32),
    String(String),
}

/// Delta encode positions
#[inline]
fn delta_encode(values: &[u64]) -> Vec<i64> {
    let mut deltas = Vec::with_capacity(values.len());
    if values.is_empty() {
        return deltas;
    }
    deltas.push(values[0] as i64);
    for i in 1..values.len() {
        deltas.push(values[i] as i64 - values[i - 1] as i64);
    }
    deltas
}

/// Delta decode back to positions
#[inline]
fn delta_decode(deltas: &[i64]) -> Vec<u64> {
    let mut values = Vec::with_capacity(deltas.len());
    if deltas.is_empty() {
        return values;
    }
    values.push(deltas[0] as u64);
    for i in 1..deltas.len() {
        values.push((values[i - 1] as i64 + deltas[i]) as u64);
    }
    values
}

impl AlignmentRecord {
    fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        write_varint(writer, self.query_name_id)?;
        write_varint(writer, self.query_start)?;
        write_varint(writer, self.query_end)?;
        writer.write_all(&[self.strand as u8])?;
        write_varint(writer, self.target_name_id)?;
        write_varint(writer, self.target_start)?;
        write_varint(writer, self.target_end)?;
        write_varint(writer, self.residue_matches)?;
        write_varint(writer, self.alignment_block_len)?;
        writer.write_all(&[self.mapping_quality])?;
        writer.write_all(&[self.tp_type.to_u8()])?;
        writer.write_all(&[complexity_metric_to_u8(&self.complexity_metric)])?;
        write_varint(writer, self.max_complexity)?;
        self.write_tracepoints(writer)?;
        write_varint(writer, self.tags.len() as u64)?;
        for tag in &self.tags {
            tag.write(writer)?;
        }
        Ok(())
    }

    fn write_raw<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        write_varint(writer, self.query_name_id)?;
        write_varint(writer, self.query_start)?;
        write_varint(writer, self.query_end)?;
        writer.write_all(&[self.strand as u8])?;
        write_varint(writer, self.target_name_id)?;
        write_varint(writer, self.target_start)?;
        write_varint(writer, self.target_end)?;
        write_varint(writer, self.residue_matches)?;
        write_varint(writer, self.alignment_block_len)?;
        writer.write_all(&[self.mapping_quality])?;
        writer.write_all(&[self.tp_type.to_u8()])?;
        writer.write_all(&[complexity_metric_to_u8(&self.complexity_metric)])?;
        write_varint(writer, self.max_complexity)?;
        self.write_tracepoints_raw(writer)?;
        write_varint(writer, self.tags.len() as u64)?;
        for tag in &self.tags {
            tag.write(writer)?;
        }
        Ok(())
    }

    fn write_tracepoints<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        match &self.tracepoints {
            TracepointData::Standard(tps) | TracepointData::Fastga(tps) => {
                write_varint(writer, tps.len() as u64)?;
                if tps.is_empty() {
                    return Ok(());
                }
                let (first_vals, second_vals): (Vec<u64>, Vec<u64>) = tps.iter().copied().unzip();

                // FastGA-aware delta encoding:
                // - FastGA: num_differences are naturally small, use raw values
                // - Standard: query_bases are incremental, use delta encoding
                let use_delta = !matches!(self.tp_type, TracepointType::Fastga);
                let first_vals_encoded = if use_delta {
                    delta_encode(&first_vals)
                } else {
                    first_vals.iter().map(|&v| v as i64).collect()
                };

                let mut first_val_buf = Vec::with_capacity(first_vals_encoded.len() * 2);
                let mut second_val_buf = Vec::with_capacity(second_vals.len() * 2);

                for &val in &first_vals_encoded {
                    let zigzag = ((val << 1) ^ (val >> 63)) as u64;
                    write_varint(&mut first_val_buf, zigzag)?;
                }
                for &val in &second_vals {
                    write_varint(&mut second_val_buf, val)?;
                }

                let first_compressed = zstd::encode_all(&first_val_buf[..], 3)?;
                let second_compressed = zstd::encode_all(&second_val_buf[..], 3)?;

                write_varint(writer, first_compressed.len() as u64)?;
                writer.write_all(&first_compressed)?;
                write_varint(writer, second_compressed.len() as u64)?;
                writer.write_all(&second_compressed)?;
            }
            TracepointData::Variable(tps) => {
                write_varint(writer, tps.len() as u64)?;
                for (a, b_opt) in tps {
                    write_varint(writer, *a)?;
                    if let Some(b) = b_opt {
                        writer.write_all(&[1])?;
                        write_varint(writer, *b)?;
                    } else {
                        writer.write_all(&[0])?;
                    }
                }
            }
            TracepointData::Mixed(items) => {
                write_varint(writer, items.len() as u64)?;
                for item in items {
                    match item {
                        MixedTracepointItem::Tracepoint(a, b) => {
                            writer.write_all(&[0])?;
                            write_varint(writer, *a)?;
                            write_varint(writer, *b)?;
                        }
                        MixedTracepointItem::CigarOp(len, op) => {
                            writer.write_all(&[1])?;
                            write_varint(writer, *len)?;
                            writer.write_all(&[*op as u8])?;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Write tracepoints WITHOUT delta encoding (for raw strategies)
    fn write_tracepoints_raw<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        match &self.tracepoints {
            TracepointData::Standard(tps) | TracepointData::Fastga(tps) => {
                write_varint(writer, tps.len() as u64)?;
                if tps.is_empty() {
                    return Ok(());
                }
                let (first_vals, second_vals): (Vec<u64>, Vec<u64>) = tps.iter().copied().unzip();

                // No delta encoding - use raw values for all types
                let mut first_val_buf = Vec::with_capacity(first_vals.len() * 2);
                let mut second_val_buf = Vec::with_capacity(second_vals.len() * 2);

                for &val in &first_vals {
                    write_varint(&mut first_val_buf, val)?;
                }
                for &val in &second_vals {
                    write_varint(&mut second_val_buf, val)?;
                }

                let first_compressed = zstd::encode_all(&first_val_buf[..], 3)?;
                let second_compressed = zstd::encode_all(&second_val_buf[..], 3)?;

                write_varint(writer, first_compressed.len() as u64)?;
                writer.write_all(&first_compressed)?;
                write_varint(writer, second_compressed.len() as u64)?;
                writer.write_all(&second_compressed)?;
            }
            TracepointData::Variable(tps) => {
                write_varint(writer, tps.len() as u64)?;
                for (a, b_opt) in tps {
                    write_varint(writer, *a)?;
                    if let Some(b) = b_opt {
                        writer.write_all(&[1])?;
                        write_varint(writer, *b)?;
                    } else {
                        writer.write_all(&[0])?;
                    }
                }
            }
            TracepointData::Mixed(items) => {
                write_varint(writer, items.len() as u64)?;
                for item in items {
                    match item {
                        MixedTracepointItem::Tracepoint(a, b) => {
                            writer.write_all(&[0])?;
                            write_varint(writer, *a)?;
                            write_varint(writer, *b)?;
                        }
                        MixedTracepointItem::CigarOp(len, op) => {
                            writer.write_all(&[1])?;
                            write_varint(writer, *len)?;
                            writer.write_all(&[*op as u8])?;
                        }
                    }
                }
            }
        }
        Ok(())
    }

    fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let query_name_id = read_varint(reader)?;
        let query_start = read_varint(reader)?;
        let query_end = read_varint(reader)?;
        let mut strand_buf = [0u8; 1];
        reader.read_exact(&mut strand_buf)?;
        let strand = strand_buf[0] as char;
        let target_name_id = read_varint(reader)?;
        let target_start = read_varint(reader)?;
        let target_end = read_varint(reader)?;
        let residue_matches = read_varint(reader)?;
        let alignment_block_len = read_varint(reader)?;
        let mut mapq_buf = [0u8; 1];
        reader.read_exact(&mut mapq_buf)?;
        let mapping_quality = mapq_buf[0];
        let mut tp_type_buf = [0u8; 1];
        reader.read_exact(&mut tp_type_buf)?;
        let tp_type = TracepointType::from_u8(tp_type_buf[0])
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let mut metric_buf = [0u8; 1];
        reader.read_exact(&mut metric_buf)?;
        let complexity_metric = complexity_metric_from_u8(metric_buf[0])?;
        let max_complexity = read_varint(reader)?;
        let tracepoints = Self::read_tracepoints(reader, tp_type)?;
        let num_tags = read_varint(reader)? as usize;
        let mut tags = Vec::with_capacity(num_tags);
        for _ in 0..num_tags {
            tags.push(Tag::read(reader)?);
        }
        Ok(Self {
            query_name_id,
            query_start,
            query_end,
            strand,
            target_name_id,
            target_start,
            target_end,
            residue_matches,
            alignment_block_len,
            mapping_quality,
            tp_type,
            complexity_metric,
            max_complexity,
            tracepoints,
            tags,
        })
    }

    fn read_tracepoints<R: Read>(
        reader: &mut R,
        tp_type: TracepointType,
    ) -> io::Result<TracepointData> {
        let num_items = read_varint(reader)? as usize;
        match tp_type {
            TracepointType::Standard | TracepointType::Fastga => {
                if num_items == 0 {
                    return Ok(match tp_type {
                        TracepointType::Standard => TracepointData::Standard(Vec::new()),
                        _ => TracepointData::Fastga(Vec::new()),
                    });
                }
                let pos_len = read_varint(reader)? as usize;
                let mut pos_compressed = vec![0u8; pos_len];
                reader.read_exact(&mut pos_compressed)?;
                let score_len = read_varint(reader)? as usize;
                let mut score_compressed = vec![0u8; score_len];
                reader.read_exact(&mut score_compressed)?;

                let pos_buf = zstd::decode_all(&pos_compressed[..])?;
                let score_buf = zstd::decode_all(&score_compressed[..])?;

                let mut pos_reader = &pos_buf[..];
                let mut pos_values = Vec::with_capacity(num_items);
                for _ in 0..num_items {
                    let zigzag = read_varint(&mut pos_reader)?;
                    let val = ((zigzag >> 1) as i64) ^ -((zigzag & 1) as i64);
                    pos_values.push(val);
                }

                // Version 1 (Varint) uses FastGA-aware delta encoding
                let positions: Vec<u64> = if matches!(tp_type, TracepointType::Fastga) {
                    // FastGA: raw values (no delta)
                    pos_values.iter().map(|&v| v as u64).collect()
                } else {
                    // Standard: delta-encoded
                    delta_decode(&pos_values)
                };

                let mut score_reader = &score_buf[..];
                let mut scores = Vec::with_capacity(num_items);
                for _ in 0..num_items {
                    scores.push(read_varint(&mut score_reader)?);
                }

                let tps: Vec<(u64, u64)> = positions.into_iter().zip(scores).collect();
                Ok(match tp_type {
                    TracepointType::Standard => TracepointData::Standard(tps),
                    _ => TracepointData::Fastga(tps),
                })
            }
            TracepointType::Variable => {
                let mut tps = Vec::with_capacity(num_items);
                for _ in 0..num_items {
                    let a = read_varint(reader)?;
                    let mut flag = [0u8; 1];
                    reader.read_exact(&mut flag)?;
                    let b_opt = if flag[0] == 1 {
                        Some(read_varint(reader)?)
                    } else {
                        None
                    };
                    tps.push((a, b_opt));
                }
                Ok(TracepointData::Variable(tps))
            }
            TracepointType::Mixed => {
                let mut items = Vec::with_capacity(num_items);
                for _ in 0..num_items {
                    let mut item_type = [0u8; 1];
                    reader.read_exact(&mut item_type)?;
                    match item_type[0] {
                        0 => {
                            let a = read_varint(reader)?;
                            let b = read_varint(reader)?;
                            items.push(MixedTracepointItem::Tracepoint(a, b));
                        }
                        1 => {
                            let len = read_varint(reader)?;
                            let mut op = [0u8; 1];
                            reader.read_exact(&mut op)?;
                            items.push(MixedTracepointItem::CigarOp(len, op[0] as char));
                        }
                        _ => {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidData,
                                "Invalid mixed item type",
                            ))
                        }
                    }
                }
                Ok(TracepointData::Mixed(items))
            }
        }
    }
}

impl Tag {
    fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.key)?;
        writer.write_all(&[self.tag_type])?;
        match &self.value {
            TagValue::Int(v) => writer.write_all(&v.to_le_bytes())?,
            TagValue::Float(v) => writer.write_all(&v.to_le_bytes())?,
            TagValue::String(s) => {
                write_varint(writer, s.len() as u64)?;
                writer.write_all(s.as_bytes())?;
            }
        }
        Ok(())
    }

    fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut key = [0u8; 2];
        reader.read_exact(&mut key)?;
        let mut tag_type = [0u8; 1];
        reader.read_exact(&mut tag_type)?;
        let value = match tag_type[0] {
            b'i' => {
                let mut buf = [0u8; 8];
                reader.read_exact(&mut buf)?;
                TagValue::Int(i64::from_le_bytes(buf))
            }
            b'f' => {
                let mut buf = [0u8; 4];
                reader.read_exact(&mut buf)?;
                TagValue::Float(f32::from_le_bytes(buf))
            }
            b'Z' => {
                let len = read_varint(reader)? as usize;
                let mut buf = vec![0u8; len];
                reader.read_exact(&mut buf)?;
                let s = String::from_utf8(buf)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                TagValue::String(s)
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid tag type",
                ))
            }
        };
        Ok(Self {
            key,
            tag_type: tag_type[0],
            value,
        })
    }
}

// ============================================================================
// PUBLIC API
// ============================================================================

/// Detect if file is binary PAF
pub fn is_binary_paf(path: &str) -> io::Result<bool> {
    if path == "-" {
        return Ok(false);
    }
    let mut file = File::open(path)?;
    let mut magic = [0u8; 4];
    match file.read_exact(&mut magic) {
        Ok(()) => Ok(&magic == BINARY_MAGIC),
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(e),
    }
}

/// Open PAF file for reading, supporting plain text, gzip, and bgzip formats
///
/// Handles three input types:
/// - `-`: Read from stdin
/// - `.gz` or `.bgz` files: Decompress with gzip decoder
/// - Plain text: Read directly
fn open_paf_reader(input_path: &str) -> io::Result<Box<dyn BufRead>> {
    if input_path == "-" {
        Ok(Box::new(BufReader::new(io::stdin())))
    } else if input_path.ends_with(".gz") || input_path.ends_with(".bgz") {
        let file = File::open(input_path).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("Failed to open input file '{}': {}", input_path, e),
            )
        })?;
        let decoder = MultiGzDecoder::new(file);
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        let file = File::open(input_path).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("Failed to open input file '{}': {}", input_path, e),
            )
        })?;
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Encode PAF with CIGAR to binary with tracepoints
/// Encode CIGAR to tracepoints and write binary format
pub fn encode_cigar_to_binary(
    input_path: &str,
    output_path: &str,
    tp_type: &TracepointType,
    max_complexity: usize,
    complexity_metric: &ComplexityMetric,
    strategy: CompressionStrategy,
) -> io::Result<()> {
    info!("Encoding CIGAR with {} strategy...", strategy);

    let input = open_paf_reader(input_path)?;
    let mut string_table = StringTable::new();
    let mut records = Vec::new();

    for (line_num, line_result) in input.lines().enumerate() {
        let line = line_result?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        match parse_paf_with_cigar(
            &line,
            &mut string_table,
            tp_type,
            max_complexity,
            complexity_metric,
        ) {
            Ok(record) => records.push(record),
            Err(e) => {
                error!("Line {}: {}", line_num + 1, e);
                return Err(e);
            }
        }
    }

    match strategy {
        CompressionStrategy::Varint => {
            write_binary(output_path, &records, &string_table)?;
        }
        CompressionStrategy::VarintRaw => {
            write_binary_raw(output_path, &records, &string_table)?;
        }
        CompressionStrategy::Huffman => {
            info!("Training Huffman codecs from {} records (FastGA-aware delta)...", records.len());
            let (first_val_codec, second_val_codec) = train_codecs_from_records(&records)?;
            info!(
                "Codecs trained: {} first_val symbols, {} second_val symbols",
                first_val_codec.num_symbols(),
                second_val_codec.num_symbols()
            );
            write_binary_adaptive(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
        }
        CompressionStrategy::HuffmanRaw => {
            info!("Training Huffman codecs from {} records (no delta)...", records.len());
            let (first_val_codec, second_val_codec) = train_codecs_from_records_raw(&records)?;
            info!(
                "Codecs trained: {} first_val symbols, {} second_val symbols",
                first_val_codec.num_symbols(),
                second_val_codec.num_symbols()
            );
            write_binary_adaptive_raw(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
        }
        CompressionStrategy::VarintHuffman => {
            info!("Training byte-level Huffman codecs from {} records (varint+huffman hybrid)...", records.len());
            let (first_val_codec, second_val_codec) = train_codecs_from_bytes(&records)?;
            info!(
                "Byte codecs trained: {} first_val byte symbols, {} second_val byte symbols",
                first_val_codec.num_symbols(),
                second_val_codec.num_symbols()
            );
            write_binary_varint_huffman(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
        }
        CompressionStrategy::VarintHuffmanRaw => {
            info!("Training byte-level Huffman codecs from {} records (varint+huffman-raw hybrid)...", records.len());
            let (first_val_codec, second_val_codec) = train_codecs_from_bytes_raw(&records)?;
            info!(
                "Byte codecs trained: {} first_val byte symbols, {} second_val byte symbols",
                first_val_codec.num_symbols(),
                second_val_codec.num_symbols()
            );
            write_binary_varint_huffman_raw(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
        }
        CompressionStrategy::Smart => {
            // Analyze data and select optimal strategy
            let result = analyze_and_select_strategy(&records)?;
            let selected_strategy = result.strategy;
            let cached_codecs = result.codecs;

            // Use selected strategy with cached codecs (if applicable)
            match selected_strategy {
                CompressionStrategy::Varint => {
                    write_binary(output_path, &records, &string_table)?;
                }
                CompressionStrategy::VarintRaw => {
                    write_binary_raw(output_path, &records, &string_table)?;
                }
                CompressionStrategy::Huffman => {
                    let (first_val_codec, second_val_codec) = cached_codecs.expect("Huffman strategy should have cached codecs");
                    info!(
                        "Using cached codecs: {} first_val symbols, {} second_val symbols",
                        first_val_codec.num_symbols(),
                        second_val_codec.num_symbols()
                    );
                    write_binary_adaptive(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
                }
                CompressionStrategy::HuffmanRaw => {
                    let (first_val_codec, second_val_codec) = cached_codecs.expect("HuffmanRaw strategy should have cached codecs");
                    info!(
                        "Using cached codecs: {} first_val symbols, {} second_val symbols",
                        first_val_codec.num_symbols(),
                        second_val_codec.num_symbols()
                    );
                    write_binary_adaptive_raw(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
                }
                CompressionStrategy::VarintHuffman => {
                    let (first_val_codec, second_val_codec) = cached_codecs.expect("VarintHuffman strategy should have cached codecs");
                    info!(
                        "Using cached byte codecs: {} first_val symbols, {} second_val symbols",
                        first_val_codec.num_symbols(),
                        second_val_codec.num_symbols()
                    );
                    write_binary_varint_huffman(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
                }
                CompressionStrategy::VarintHuffmanRaw => {
                    let (first_val_codec, second_val_codec) = cached_codecs.expect("VarintHuffmanRaw strategy should have cached codecs");
                    info!(
                        "Using cached byte codecs: {} first_val symbols, {} second_val symbols",
                        first_val_codec.num_symbols(),
                        second_val_codec.num_symbols()
                    );
                    write_binary_varint_huffman_raw(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
                }
                CompressionStrategy::Smart => {
                    // Should never happen, but prevent infinite recursion
                    unreachable!("Smart strategy should not select itself")
                }
            }
        }
    }

    info!(
        "Encoded {} records ({} unique names) with {} strategy",
        records.len(),
        string_table.len(),
        strategy
    );
    Ok(())
}

/// Convert binary PAF to text format
pub fn decompress_paf(input_path: &str, output_path: &str) -> io::Result<()> {
    info!("Decompressing {} to text format...", input_path);

    let input = File::open(input_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to open input file '{}': {}", input_path, e),
        )
    })?;
    let mut reader = BufReader::new(input);

    let header = BinaryPafHeader::read(&mut reader)?;

    if header.version != 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Unsupported format version: {}", header.version),
        ));
    }

    let strategy = CompressionStrategy::from_code(header.flags)?;
    info!(
        "Reading {} records ({} unique names) [{}]",
        header.num_records, header.num_strings, strategy
    );

    match strategy {
        CompressionStrategy::Varint => decompress_varint(reader, output_path, &header),
        CompressionStrategy::VarintRaw => decompress_varint_raw(reader, output_path, &header),
        CompressionStrategy::Huffman => decompress_huffman(reader, output_path, &header),
        CompressionStrategy::HuffmanRaw => decompress_huffman_raw(reader, output_path, &header),
        CompressionStrategy::VarintHuffman => decompress_varint_huffman_hybrid(reader, output_path, &header),
        CompressionStrategy::VarintHuffmanRaw => decompress_varint_huffman_raw_hybrid(reader, output_path, &header),
        CompressionStrategy::Smart => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Smart strategy should not be stored in file",
        )),
    }
}

// Version 1: Varint with FastGA-aware delta encoding
fn decompress_varint<R: Read>(
    mut reader: R,
    output_path: &str,
    header: &BinaryPafHeader,
) -> io::Result<()> {
    let mut records = Vec::with_capacity(header.num_records as usize);
    for _ in 0..header.num_records {
        records.push(AlignmentRecord::read(&mut reader)?);
    }
    let string_table = StringTable::read(&mut reader)?;
    write_paf_output(output_path, &records, &string_table)?;
    info!("Decompressed {} records", header.num_records);
    Ok(())
}

// Version 2: Huffman with FastGA-aware delta encoding
fn decompress_huffman<R: Read>(
    mut reader: R,
    output_path: &str,
    header: &BinaryPafHeader,
) -> io::Result<()> {
    let first_val_codec = AdaptiveCodec::read(&mut reader)?;
    let second_val_codec = AdaptiveCodec::read(&mut reader)?;
    info!(
        "Codecs: {} pos symbols, {} score symbols",
        first_val_codec.num_symbols(),
        second_val_codec.num_symbols()
    );

    let mut records = Vec::with_capacity(header.num_records as usize);
    for _ in 0..header.num_records {
        records.push(read_record_huffman(&mut reader, &first_val_codec, &second_val_codec, true)?);
    }
    let string_table = StringTable::read(&mut reader)?;
    write_paf_output(output_path, &records, &string_table)?;
    info!("Decompressed {} records", header.num_records);
    Ok(())
}

// Version 3: Huffman without delta encoding
fn decompress_huffman_raw<R: Read>(
    mut reader: R,
    output_path: &str,
    header: &BinaryPafHeader,
) -> io::Result<()> {
    let first_val_codec = AdaptiveCodec::read(&mut reader)?;
    let second_val_codec = AdaptiveCodec::read(&mut reader)?;
    let mut records = Vec::with_capacity(header.num_records as usize);
    for _ in 0..header.num_records {
        records.push(read_record_huffman(&mut reader, &first_val_codec, &second_val_codec, false)?);
    }
    let string_table = StringTable::read(&mut reader)?;
    write_paf_output(output_path, &records, &string_table)?;
    info!("Decompressed {} records", header.num_records);
    Ok(())
}

// Version 4: Varint without delta encoding
fn decompress_varint_raw<R: Read>(
    mut reader: R,
    output_path: &str,
    header: &BinaryPafHeader,
) -> io::Result<()> {
    let mut records = Vec::with_capacity(header.num_records as usize);
    for _ in 0..header.num_records {
        records.push(read_record_varint(&mut reader, false)?);
    }
    let string_table = StringTable::read(&mut reader)?;
    write_paf_output(output_path, &records, &string_table)?;
    info!("Decompressed {} records", header.num_records);
    Ok(())
}

// Version 5: Varint+Huffman hybrid with delta
fn decompress_varint_huffman_hybrid<R: Read>(
    mut reader: R,
    output_path: &str,
    header: &BinaryPafHeader,
) -> io::Result<()> {
    let first_val_codec = AdaptiveCodec::read(&mut reader)?;
    let second_val_codec = AdaptiveCodec::read(&mut reader)?;
    let mut records = Vec::with_capacity(header.num_records as usize);
    for _ in 0..header.num_records {
        records.push(read_record_varint_huffman(&mut reader, &first_val_codec, &second_val_codec, true)?);
    }
    let string_table = StringTable::read(&mut reader)?;
    write_paf_output(output_path, &records, &string_table)?;
    info!("Decompressed {} records", header.num_records);
    Ok(())
}

// Version 6: Varint+Huffman hybrid without delta
fn decompress_varint_huffman_raw_hybrid<R: Read>(
    mut reader: R,
    output_path: &str,
    header: &BinaryPafHeader,
) -> io::Result<()> {
    let first_val_codec = AdaptiveCodec::read(&mut reader)?;
    let second_val_codec = AdaptiveCodec::read(&mut reader)?;
    let mut records = Vec::with_capacity(header.num_records as usize);
    for _ in 0..header.num_records {
        records.push(read_record_varint_huffman(&mut reader, &first_val_codec, &second_val_codec, false)?);
    }
    let string_table = StringTable::read(&mut reader)?;
    write_paf_output(output_path, &records, &string_table)?;
    info!("Decompressed {} records", header.num_records);
    Ok(())
}

fn read_record_huffman<R: Read>(
    reader: &mut R,
    first_val_codec: &AdaptiveCodec,
    second_val_codec: &AdaptiveCodec,
    fastga_aware_delta: bool,
) -> io::Result<AlignmentRecord> {
    let query_name_id = read_varint(reader)?;
    let query_start = read_varint(reader)?;
    let query_end = read_varint(reader)?;
    let mut strand_buf = [0u8; 1];
    reader.read_exact(&mut strand_buf)?;
    let strand = strand_buf[0] as char;
    let target_name_id = read_varint(reader)?;
    let target_start = read_varint(reader)?;
    let target_end = read_varint(reader)?;
    let residue_matches = read_varint(reader)?;
    let alignment_block_len = read_varint(reader)?;
    let mut mapq_buf = [0u8; 1];
    reader.read_exact(&mut mapq_buf)?;
    let mapping_quality = mapq_buf[0];
    let mut tp_type_buf = [0u8; 1];
    reader.read_exact(&mut tp_type_buf)?;
    let tp_type = TracepointType::from_u8(tp_type_buf[0])
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut metric_buf = [0u8; 1];
    reader.read_exact(&mut metric_buf)?;
    let complexity_metric = complexity_metric_from_u8(metric_buf[0])?;
    let max_complexity = read_varint(reader)?;

    let num_items = read_varint(reader)? as usize;
    let tracepoints = if num_items == 0 {
        match tp_type {
            TracepointType::Standard => TracepointData::Standard(Vec::new()),
            _ => TracepointData::Fastga(Vec::new()),
        }
    } else {
        let pos_len = read_varint(reader)? as usize;
        let mut pos_compressed = vec![0u8; pos_len];
        reader.read_exact(&mut pos_compressed)?;
        let score_len = read_varint(reader)? as usize;
        let mut score_compressed = vec![0u8; score_len];
        reader.read_exact(&mut score_compressed)?;

        let pos_bytes = zstd::decode_all(&pos_compressed[..])?;
        let score_bytes = zstd::decode_all(&score_compressed[..])?;

        let mut pos_reader = BitReader::new(&pos_bytes);
        let mut pos_values = Vec::with_capacity(num_items);
        for _ in 0..num_items {
            pos_values.push(first_val_codec.decode(&mut pos_reader)?);
        }

        // Apply delta decoding based on strategy and tracepoint type
        let positions: Vec<u64> = if fastga_aware_delta {
            // FastGA-aware: only delta-decode for Standard, use raw for FastGA
            if matches!(tp_type, TracepointType::Fastga) {
                pos_values.iter().map(|&v| v as u64).collect()
            } else {
                delta_decode(&pos_values)
            }
        } else {
            // Raw strategy: never delta-decode
            pos_values.iter().map(|&v| v as u64).collect()
        };

        let mut score_reader = BitReader::new(&score_bytes);
        let mut scores = Vec::with_capacity(num_items);
        for _ in 0..num_items {
            scores.push(second_val_codec.decode(&mut score_reader)? as u64);
        }

        let tps: Vec<(u64, u64)> = positions.into_iter().zip(scores).collect();
        match tp_type {
            TracepointType::Standard => TracepointData::Standard(tps),
            _ => TracepointData::Fastga(tps),
        }
    };

    let num_tags = read_varint(reader)? as usize;
    let mut tags = Vec::with_capacity(num_tags);
    for _ in 0..num_tags {
        tags.push(Tag::read(reader)?);
    }

    Ok(AlignmentRecord {
        query_name_id,
        query_start,
        query_end,
        strand,
        target_name_id,
        target_start,
        target_end,
        residue_matches,
        alignment_block_len,
        mapping_quality,
        tp_type,
        complexity_metric,
        max_complexity,
        tracepoints,
        tags,
    })
}

// Read record for version 4 (VarintRaw) - always use raw values, no delta
fn read_record_varint<R: Read>(
    reader: &mut R,
    _use_delta: bool, // Parameter for consistency, always false for VarintRaw
) -> io::Result<AlignmentRecord> {
    // Read PAF fields (same as version 1)
    let query_name_id = read_varint(reader)?;
    let query_start = read_varint(reader)?;
    let query_end = read_varint(reader)?;
    let mut strand_buf = [0u8; 1];
    reader.read_exact(&mut strand_buf)?;
    let strand = strand_buf[0] as char;
    let target_name_id = read_varint(reader)?;
    let target_start = read_varint(reader)?;
    let target_end = read_varint(reader)?;
    let residue_matches = read_varint(reader)?;
    let alignment_block_len = read_varint(reader)?;
    let mut mapq_buf = [0u8; 1];
    reader.read_exact(&mut mapq_buf)?;
    let mapping_quality = mapq_buf[0];
    let mut tp_type_buf = [0u8; 1];
    reader.read_exact(&mut tp_type_buf)?;
    let tp_type = TracepointType::from_u8(tp_type_buf[0])
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut metric_buf = [0u8; 1];
    reader.read_exact(&mut metric_buf)?;
    let complexity_metric = complexity_metric_from_u8(metric_buf[0])?;
    let max_complexity = read_varint(reader)?;

    // Read tracepoints (raw, no delta)
    let tracepoints = read_tracepoints_raw(reader, tp_type)?;

    let num_tags = read_varint(reader)? as usize;
    let mut tags = Vec::with_capacity(num_tags);
    for _ in 0..num_tags {
        tags.push(Tag::read(reader)?);
    }

    Ok(AlignmentRecord {
        query_name_id,
        query_start,
        query_end,
        strand,
        target_name_id,
        target_start,
        target_end,
        residue_matches,
        alignment_block_len,
        mapping_quality,
        tp_type,
        complexity_metric,
        max_complexity,
        tracepoints,
        tags,
    })
}

// Read tracepoints without delta encoding (version 4)
fn read_tracepoints_raw<R: Read>(
    reader: &mut R,
    tp_type: TracepointType,
) -> io::Result<TracepointData> {
    let num_items = read_varint(reader)? as usize;
    if num_items == 0 {
        return Ok(match tp_type {
            TracepointType::Standard => TracepointData::Standard(Vec::new()),
            _ => TracepointData::Fastga(Vec::new()),
        });
    }

    let pos_len = read_varint(reader)? as usize;
    let mut pos_compressed = vec![0u8; pos_len];
    reader.read_exact(&mut pos_compressed)?;
    let score_len = read_varint(reader)? as usize;
    let mut score_compressed = vec![0u8; score_len];
    reader.read_exact(&mut score_compressed)?;

    let pos_buf = zstd::decode_all(&pos_compressed[..])?;
    let score_buf = zstd::decode_all(&score_compressed[..])?;

    // Read raw varint values (no zigzag, no delta)
    let mut pos_reader = &pos_buf[..];
    let mut positions = Vec::with_capacity(num_items);
    for _ in 0..num_items {
        positions.push(read_varint(&mut pos_reader)?);
    }

    let mut score_reader = &score_buf[..];
    let mut scores = Vec::with_capacity(num_items);
    for _ in 0..num_items {
        scores.push(read_varint(&mut score_reader)?);
    }

    let tps: Vec<(u64, u64)> = positions.into_iter().zip(scores).collect();
    Ok(match tp_type {
        TracepointType::Standard => TracepointData::Standard(tps),
        _ => TracepointData::Fastga(tps),
    })
}

// Read record for versions 5 & 6 (VarintHuffman variants)
fn read_record_varint_huffman<R: Read>(
    reader: &mut R,
    first_val_codec: &AdaptiveCodec,
    second_val_codec: &AdaptiveCodec,
    fastga_aware_delta: bool,
) -> io::Result<AlignmentRecord> {
    // Same PAF field reading as Huffman
    let query_name_id = read_varint(reader)?;
    let query_start = read_varint(reader)?;
    let query_end = read_varint(reader)?;
    let mut strand_buf = [0u8; 1];
    reader.read_exact(&mut strand_buf)?;
    let strand = strand_buf[0] as char;
    let target_name_id = read_varint(reader)?;
    let target_start = read_varint(reader)?;
    let target_end = read_varint(reader)?;
    let residue_matches = read_varint(reader)?;
    let alignment_block_len = read_varint(reader)?;
    let mut mapq_buf = [0u8; 1];
    reader.read_exact(&mut mapq_buf)?;
    let mapping_quality = mapq_buf[0];
    let mut tp_type_buf = [0u8; 1];
    reader.read_exact(&mut tp_type_buf)?;
    let tp_type = TracepointType::from_u8(tp_type_buf[0])
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let mut metric_buf = [0u8; 1];
    reader.read_exact(&mut metric_buf)?;
    let complexity_metric = complexity_metric_from_u8(metric_buf[0])?;
    let max_complexity = read_varint(reader)?;

    let num_items = read_varint(reader)? as usize;
    let tracepoints = if num_items == 0 {
        match tp_type {
            TracepointType::Standard => TracepointData::Standard(Vec::new()),
            _ => TracepointData::Fastga(Vec::new()),
        }
    } else {
        // Varint+Huffman: varint → bytes → Huffman → zstd
        let pos_len = read_varint(reader)? as usize;
        let mut pos_compressed = vec![0u8; pos_len];
        reader.read_exact(&mut pos_compressed)?;
        let score_len = read_varint(reader)? as usize;
        let mut score_compressed = vec![0u8; score_len];
        reader.read_exact(&mut score_compressed)?;

        let pos_bytes = zstd::decode_all(&pos_compressed[..])?;
        let score_bytes = zstd::decode_all(&score_compressed[..])?;

        // Decode bytes with Huffman - decode until no more data
        let mut pos_reader = BitReader::new(&pos_bytes);
        let mut pos_byte_values = Vec::new();
        loop {
            match first_val_codec.decode(&mut pos_reader) {
                Ok(byte_val) => pos_byte_values.push(byte_val as u8),
                Err(_) => break,
            }
        }

        // Decode varint from bytes (zigzag for v5, plain for v6)
        let mut byte_reader = &pos_byte_values[..];
        let mut pos_values = Vec::with_capacity(num_items);
        for _ in 0..num_items {
            if fastga_aware_delta {
                // Version 5: zigzag encoding
                let zigzag = read_varint(&mut byte_reader)?;
                let val = ((zigzag >> 1) as i64) ^ -((zigzag & 1) as i64);
                pos_values.push(val);
            } else {
                // Version 6: plain varint
                pos_values.push(read_varint(&mut byte_reader)? as i64);
            }
        }

        // Apply delta based on strategy
        let positions: Vec<u64> = if fastga_aware_delta {
            if matches!(tp_type, TracepointType::Fastga) {
                pos_values.iter().map(|&v| v as u64).collect()
            } else {
                delta_decode(&pos_values)
            }
        } else {
            pos_values.iter().map(|&v| v as u64).collect()
        };

        // Decode second values - decode until no more data
        let mut score_reader = BitReader::new(&score_bytes);
        let mut score_byte_values = Vec::new();
        loop {
            match second_val_codec.decode(&mut score_reader) {
                Ok(byte_val) => score_byte_values.push(byte_val as u8),
                Err(_) => break,
            }
        }

        let mut byte_reader = &score_byte_values[..];
        let mut scores = Vec::with_capacity(num_items);
        for _ in 0..num_items {
            scores.push(read_varint(&mut byte_reader)?);
        }

        let tps: Vec<(u64, u64)> = positions.into_iter().zip(scores).collect();
        match tp_type {
            TracepointType::Standard => TracepointData::Standard(tps),
            _ => TracepointData::Fastga(tps),
        }
    };

    let num_tags = read_varint(reader)? as usize;
    let mut tags = Vec::with_capacity(num_tags);
    for _ in 0..num_tags {
        tags.push(Tag::read(reader)?);
    }

    Ok(AlignmentRecord {
        query_name_id,
        query_start,
        query_end,
        strand,
        target_name_id,
        target_start,
        target_end,
        residue_matches,
        alignment_block_len,
        mapping_quality,
        tp_type,
        complexity_metric,
        max_complexity,
        tracepoints,
        tags,
    })
}

/// Compress PAF with tracepoints to binary format
///
/// Uses delta encoding + varint + zstd compression for optimal balance
/// of speed and compression ratio on genomic alignment data.
/// Compress PAF with tracepoints to binary format
pub fn compress_paf(input_path: &str, output_path: &str, strategy: CompressionStrategy) -> io::Result<()> {
    info!("Compressing PAF with {} strategy...", strategy);

    let input = open_paf_reader(input_path)?;
    let mut string_table = StringTable::new();
    let mut records = Vec::new();

    for (line_num, line_result) in input.lines().enumerate() {
        let line = line_result?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        match parse_paf_with_tracepoints(&line, &mut string_table) {
            Ok(record) => records.push(record),
            Err(e) => {
                error!("Line {}: {}", line_num + 1, e);
                return Err(e);
            }
        }
    }

    match strategy {
        CompressionStrategy::Varint => {
            write_binary(output_path, &records, &string_table)?;
        }
        CompressionStrategy::VarintRaw => {
            write_binary_raw(output_path, &records, &string_table)?;
        }
        CompressionStrategy::Huffman => {
            info!("Training Huffman codecs from {} records (FastGA-aware delta)...", records.len());
            let (first_val_codec, second_val_codec) = train_codecs_from_records(&records)?;
            info!(
                "Codecs trained: {} first_val symbols, {} second_val symbols",
                first_val_codec.num_symbols(),
                second_val_codec.num_symbols()
            );
            write_binary_adaptive(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
        }
        CompressionStrategy::HuffmanRaw => {
            info!("Training Huffman codecs from {} records (no delta)...", records.len());
            let (first_val_codec, second_val_codec) = train_codecs_from_records_raw(&records)?;
            info!(
                "Codecs trained: {} first_val symbols, {} second_val symbols",
                first_val_codec.num_symbols(),
                second_val_codec.num_symbols()
            );
            write_binary_adaptive_raw(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
        }
        CompressionStrategy::VarintHuffman => {
            info!("Training byte-level Huffman codecs from {} records (varint+huffman hybrid)...", records.len());
            let (first_val_codec, second_val_codec) = train_codecs_from_bytes(&records)?;
            info!(
                "Byte codecs trained: {} first_val byte symbols, {} second_val byte symbols",
                first_val_codec.num_symbols(),
                second_val_codec.num_symbols()
            );
            write_binary_varint_huffman(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
        }
        CompressionStrategy::VarintHuffmanRaw => {
            info!("Training byte-level Huffman codecs from {} records (varint+huffman-raw hybrid)...", records.len());
            let (first_val_codec, second_val_codec) = train_codecs_from_bytes_raw(&records)?;
            info!(
                "Byte codecs trained: {} first_val byte symbols, {} second_val byte symbols",
                first_val_codec.num_symbols(),
                second_val_codec.num_symbols()
            );
            write_binary_varint_huffman_raw(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
        }
        CompressionStrategy::Smart => {
            // Analyze data and select optimal strategy
            let result = analyze_and_select_strategy(&records)?;
            let selected_strategy = result.strategy;
            let cached_codecs = result.codecs;

            // Use selected strategy with cached codecs (if applicable)
            match selected_strategy {
                CompressionStrategy::Varint => {
                    write_binary(output_path, &records, &string_table)?;
                }
                CompressionStrategy::VarintRaw => {
                    write_binary_raw(output_path, &records, &string_table)?;
                }
                CompressionStrategy::Huffman => {
                    let (first_val_codec, second_val_codec) = cached_codecs.expect("Huffman strategy should have cached codecs");
                    info!(
                        "Using cached codecs: {} first_val symbols, {} second_val symbols",
                        first_val_codec.num_symbols(),
                        second_val_codec.num_symbols()
                    );
                    write_binary_adaptive(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
                }
                CompressionStrategy::HuffmanRaw => {
                    let (first_val_codec, second_val_codec) = cached_codecs.expect("HuffmanRaw strategy should have cached codecs");
                    info!(
                        "Using cached codecs: {} first_val symbols, {} second_val symbols",
                        first_val_codec.num_symbols(),
                        second_val_codec.num_symbols()
                    );
                    write_binary_adaptive_raw(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
                }
                CompressionStrategy::VarintHuffman => {
                    let (first_val_codec, second_val_codec) = cached_codecs.expect("VarintHuffman strategy should have cached codecs");
                    info!(
                        "Using cached byte codecs: {} first_val symbols, {} second_val symbols",
                        first_val_codec.num_symbols(),
                        second_val_codec.num_symbols()
                    );
                    write_binary_varint_huffman(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
                }
                CompressionStrategy::VarintHuffmanRaw => {
                    let (first_val_codec, second_val_codec) = cached_codecs.expect("VarintHuffmanRaw strategy should have cached codecs");
                    info!(
                        "Using cached byte codecs: {} first_val symbols, {} second_val symbols",
                        first_val_codec.num_symbols(),
                        second_val_codec.num_symbols()
                    );
                    write_binary_varint_huffman_raw(output_path, &records, &string_table, &first_val_codec, &second_val_codec)?;
                }
                CompressionStrategy::Smart => {
                    // Should never happen, but prevent infinite recursion
                    unreachable!("Smart strategy should not select itself")
                }
            }
        }
    }

    info!(
        "Compressed {} records ({} unique names) with {} strategy",
        records.len(),
        string_table.len(),
        strategy
    );
    Ok(())
}

/// Train Huffman codecs from alignment records
///
/// Tracepoint pairs are:
/// - Standard/Mixed/Variable: (query_bases, target_bases)
/// - FastGA: (num_differences, target_bases)
///
/// Compression strategy:
/// - First value: delta encoded (difference from previous first value)
/// - Second value: raw value (target_bases in all types)
fn train_codecs_from_records(
    records: &[AlignmentRecord],
) -> io::Result<(AdaptiveCodec, AdaptiveCodec)> {
    let mut first_val_frequencies: HashMap<i64, u64> = HashMap::new();
    let mut second_val_frequencies: HashMap<i64, u64> = HashMap::new();

    // Collect frequencies from all records
    for record in records {
        let tracepoints = match &record.tracepoints {
            TracepointData::Standard(tps) | TracepointData::Fastga(tps) => tps,
            _ => continue,
        };

        if tracepoints.is_empty() {
            continue;
        }

        // FastGA-aware delta encoding:
        // - FastGA: num_differences are naturally small, use raw values
        // - Standard: query_bases are incremental, use delta encoding
        let use_delta = !matches!(record.tp_type, TracepointType::Fastga);

        if use_delta {
            // Delta encoding for Standard type
            let mut prev_first = 0i64;
            for &(first_val, second_val) in tracepoints {
                let delta = first_val as i64 - prev_first;
                *first_val_frequencies.entry(delta).or_insert(0) += 1;
                prev_first = first_val as i64;

                *second_val_frequencies.entry(second_val as i64).or_insert(0) += 1;
            }
        } else {
            // Raw values for FastGA type
            for &(first_val, second_val) in tracepoints {
                *first_val_frequencies.entry(first_val as i64).or_insert(0) += 1;
                *second_val_frequencies.entry(second_val as i64).or_insert(0) += 1;
            }
        }
    }

    if first_val_frequencies.is_empty() || second_val_frequencies.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "No tracepoint data found for training codecs",
        ));
    }

    // Debug: log frequency distribution statistics
    debug!("First value delta statistics:");
    debug!("  Unique symbols: {}", first_val_frequencies.len());
    let first_total: u64 = first_val_frequencies.values().sum();
    let mut first_sorted: Vec<_> = first_val_frequencies.iter().collect();
    first_sorted.sort_by(|a, b| b.1.cmp(a.1));
    if let Some((val, freq)) = first_sorted.first() {
        debug!("  Most frequent: {} ({}x, {:.2}%)", val, freq, (**freq as f64 / first_total as f64) * 100.0);
    }
    if first_sorted.len() >= 3 {
        let top3_freq: u64 = first_sorted.iter().take(3).map(|(_, f)| **f).sum();
        debug!("  Top 3 symbols: {:.2}% of data", (top3_freq as f64 / first_total as f64) * 100.0);
    }

    debug!("Second value statistics:");
    debug!("  Unique symbols: {}", second_val_frequencies.len());
    let second_total: u64 = second_val_frequencies.values().sum();
    let mut second_sorted: Vec<_> = second_val_frequencies.iter().collect();
    second_sorted.sort_by(|a, b| b.1.cmp(a.1));
    if let Some((val, freq)) = second_sorted.first() {
        debug!("  Most frequent: {} ({}x, {:.2}%)", val, freq, (**freq as f64 / second_total as f64) * 100.0);
    }
    if second_sorted.len() >= 3 {
        let top3_freq: u64 = second_sorted.iter().take(3).map(|(_, f)| **f).sum();
        debug!("  Top 3 symbols: {:.2}% of data", (top3_freq as f64 / second_total as f64) * 100.0);
    }

    let first_val_codec = AdaptiveCodec::build(&first_val_frequencies)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let second_val_codec = AdaptiveCodec::build(&second_val_frequencies)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    // Debug: calculate bit usage efficiency
    let first_theoretical_bits = (first_val_frequencies.len() as f64).log2().ceil();
    let mut first_weighted_bits = 0.0;
    if let Some(ref tree) = first_val_codec.huffman_tree {
        for (val, freq) in &first_val_frequencies {
            if let Some(code) = tree.encode_table.get(val) {
                first_weighted_bits += code.bits.len() as f64 * (*freq as f64);
            }
        }
    }
    let first_avg_bits = first_weighted_bits / first_total as f64;
    debug!("  Theoretical bits/symbol: {:.2}", first_theoretical_bits);
    debug!("  Huffman bits/symbol: {:.2}", first_avg_bits);
    debug!("  Bit efficiency: {:.1}%", (first_avg_bits / first_theoretical_bits) * 100.0);

    let second_theoretical_bits = (second_val_frequencies.len() as f64).log2().ceil();
    let mut second_weighted_bits = 0.0;
    if let Some(ref tree) = second_val_codec.huffman_tree {
        for (val, freq) in &second_val_frequencies {
            if let Some(code) = tree.encode_table.get(val) {
                second_weighted_bits += code.bits.len() as f64 * (*freq as f64);
            }
        }
    }
    let second_avg_bits = second_weighted_bits / second_total as f64;
    debug!("  Theoretical bits/symbol: {:.2}", second_theoretical_bits);
    debug!("  Huffman bits/symbol: {:.2}", second_avg_bits);
    debug!("  Bit efficiency: {:.1}%", (second_avg_bits / second_theoretical_bits) * 100.0);

    Ok((first_val_codec, second_val_codec))
}

/// Train Huffman codecs from alignment records (NO delta encoding)
///
/// Compression strategy:
/// - First value: raw value (no delta)
/// - Second value: raw value (target_bases in all types)
fn train_codecs_from_records_raw(
    records: &[AlignmentRecord],
) -> io::Result<(AdaptiveCodec, AdaptiveCodec)> {
    let mut first_val_frequencies: HashMap<i64, u64> = HashMap::new();
    let mut second_val_frequencies: HashMap<i64, u64> = HashMap::new();

    // Collect frequencies from all records
    for record in records {
        let tracepoints = match &record.tracepoints {
            TracepointData::Standard(tps) | TracepointData::Fastga(tps) => tps,
            _ => continue,
        };

        if tracepoints.is_empty() {
            continue;
        }

        // No delta encoding - use raw values for all types
        for &(first_val, second_val) in tracepoints {
            *first_val_frequencies.entry(first_val as i64).or_insert(0) += 1;
            *second_val_frequencies.entry(second_val as i64).or_insert(0) += 1;
        }
    }

    if first_val_frequencies.is_empty() || second_val_frequencies.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "No tracepoint data found for training codecs",
        ));
    }

    // Debug: log frequency distribution statistics
    debug!("First value (raw) statistics:");
    debug!("  Unique symbols: {}", first_val_frequencies.len());
    let first_total: u64 = first_val_frequencies.values().sum();
    let mut first_sorted: Vec<_> = first_val_frequencies.iter().collect();
    first_sorted.sort_by(|a, b| b.1.cmp(a.1));
    if let Some((val, freq)) = first_sorted.first() {
        debug!("  Most frequent: {} ({}x, {:.2}%)", val, freq, (**freq as f64 / first_total as f64) * 100.0);
    }
    if first_sorted.len() >= 3 {
        let top3_freq: u64 = first_sorted.iter().take(3).map(|(_, f)| **f).sum();
        debug!("  Top 3 symbols: {:.2}% of data", (top3_freq as f64 / first_total as f64) * 100.0);
    }

    debug!("Second value statistics:");
    debug!("  Unique symbols: {}", second_val_frequencies.len());
    let second_total: u64 = second_val_frequencies.values().sum();
    let mut second_sorted: Vec<_> = second_val_frequencies.iter().collect();
    second_sorted.sort_by(|a, b| b.1.cmp(a.1));
    if let Some((val, freq)) = second_sorted.first() {
        debug!("  Most frequent: {} ({}x, {:.2}%)", val, freq, (**freq as f64 / second_total as f64) * 100.0);
    }
    if second_sorted.len() >= 3 {
        let top3_freq: u64 = second_sorted.iter().take(3).map(|(_, f)| **f).sum();
        debug!("  Top 3 symbols: {:.2}% of data", (top3_freq as f64 / second_total as f64) * 100.0);
    }

    let first_val_codec = AdaptiveCodec::build(&first_val_frequencies)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let second_val_codec = AdaptiveCodec::build(&second_val_frequencies)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    // Debug: calculate bit usage efficiency
    let first_theoretical_bits = (first_val_frequencies.len() as f64).log2().ceil();
    let mut first_weighted_bits = 0.0;
    if let Some(ref tree) = first_val_codec.huffman_tree {
        for (val, freq) in &first_val_frequencies {
            if let Some(code) = tree.encode_table.get(val) {
                first_weighted_bits += code.bits.len() as f64 * (*freq as f64);
            }
        }
    }
    let first_avg_bits = first_weighted_bits / first_total as f64;
    debug!("  Theoretical bits/symbol: {:.2}", first_theoretical_bits);
    debug!("  Huffman bits/symbol: {:.2}", first_avg_bits);
    debug!("  Bit efficiency: {:.1}%", (first_avg_bits / first_theoretical_bits) * 100.0);

    let second_theoretical_bits = (second_val_frequencies.len() as f64).log2().ceil();
    let mut second_weighted_bits = 0.0;
    if let Some(ref tree) = second_val_codec.huffman_tree {
        for (val, freq) in &second_val_frequencies {
            if let Some(code) = tree.encode_table.get(val) {
                second_weighted_bits += code.bits.len() as f64 * (*freq as f64);
            }
        }
    }
    let second_avg_bits = second_weighted_bits / second_total as f64;
    debug!("  Theoretical bits/symbol: {:.2}", second_theoretical_bits);
    debug!("  Huffman bits/symbol: {:.2}", second_avg_bits);
    debug!("  Bit efficiency: {:.1}%", (second_avg_bits / second_theoretical_bits) * 100.0);

    Ok((first_val_codec, second_val_codec))
}

/// Train byte-level Huffman codecs by varint encoding first, then Huffman on bytes
///
/// Pipeline: Delta/Raw → Varint → Byte frequencies → Byte-level Huffman
/// Compression strategy:
/// - First value: FastGA-aware (delta for Standard, raw for FastGA)
/// - Second value: raw target_bases
/// - Huffman trained on byte distribution (0-255) from varint output
fn train_codecs_from_bytes(
    records: &[AlignmentRecord],
) -> io::Result<(AdaptiveCodec, AdaptiveCodec)> {
    let mut first_val_byte_frequencies: HashMap<i64, u64> = HashMap::new();
    let mut second_val_byte_frequencies: HashMap<i64, u64> = HashMap::new();

    // Collect byte frequencies by varint encoding all values
    for record in records {
        let tracepoints = match &record.tracepoints {
            TracepointData::Standard(tps) | TracepointData::Fastga(tps) => tps,
            _ => continue,
        };

        if tracepoints.is_empty() {
            continue;
        }

        // FastGA-aware delta encoding
        let use_delta = !matches!(record.tp_type, TracepointType::Fastga);

        // Encode first values to bytes
        let first_vals: Vec<u64> = tracepoints.iter().map(|&(v, _)| v).collect();
        let first_vals_encoded: Vec<i64> = if use_delta {
            delta_encode(&first_vals)
        } else {
            first_vals.iter().map(|&v| v as i64).collect()
        };

        // Varint encode and collect byte frequencies
        for val in first_vals_encoded {
            let mut buf = Vec::new();
            let zigzag = ((val << 1) ^ (val >> 63)) as u64;
            write_varint(&mut buf, zigzag).unwrap();
            for byte in buf {
                *first_val_byte_frequencies.entry(byte as i64).or_insert(0) += 1;
            }
        }

        // Encode second values (always raw)
        for &(_, val) in tracepoints {
            let mut buf = Vec::new();
            write_varint(&mut buf, val).unwrap();
            for byte in buf {
                *second_val_byte_frequencies.entry(byte as i64).or_insert(0) += 1;
            }
        }
    }

    if first_val_byte_frequencies.is_empty() || second_val_byte_frequencies.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "No tracepoint data found for training byte codecs",
        ));
    }

    // Debug: log byte frequency distribution statistics
    debug!("First value byte statistics (after varint):");
    debug!("  Unique byte symbols: {}", first_val_byte_frequencies.len());
    let first_total: u64 = first_val_byte_frequencies.values().sum();
    let mut first_sorted: Vec<_> = first_val_byte_frequencies.iter().collect();
    first_sorted.sort_by(|a, b| b.1.cmp(a.1));
    if let Some((val, freq)) = first_sorted.first() {
        debug!("  Most frequent: 0x{:02X} ({}x, {:.2}%)", val, freq, (**freq as f64 / first_total as f64) * 100.0);
    }
    if first_sorted.len() >= 3 {
        let top3_freq: u64 = first_sorted.iter().take(3).map(|(_, f)| **f).sum();
        debug!("  Top 3 bytes: {:.2}% of data", (top3_freq as f64 / first_total as f64) * 100.0);
    }

    debug!("Second value byte statistics (after varint):");
    debug!("  Unique byte symbols: {}", second_val_byte_frequencies.len());
    let second_total: u64 = second_val_byte_frequencies.values().sum();
    let mut second_sorted: Vec<_> = second_val_byte_frequencies.iter().collect();
    second_sorted.sort_by(|a, b| b.1.cmp(a.1));
    if let Some((val, freq)) = second_sorted.first() {
        debug!("  Most frequent: 0x{:02X} ({}x, {:.2}%)", val, freq, (**freq as f64 / second_total as f64) * 100.0);
    }
    if second_sorted.len() >= 3 {
        let top3_freq: u64 = second_sorted.iter().take(3).map(|(_, f)| **f).sum();
        debug!("  Top 3 bytes: {:.2}% of data", (top3_freq as f64 / second_total as f64) * 100.0);
    }

    let first_val_codec = AdaptiveCodec::build(&first_val_byte_frequencies)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let second_val_codec = AdaptiveCodec::build(&second_val_byte_frequencies)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((first_val_codec, second_val_codec))
}

/// Train byte-level Huffman codecs (NO delta encoding)
///
/// Pipeline: Raw → Varint → Byte frequencies → Byte-level Huffman
fn train_codecs_from_bytes_raw(
    records: &[AlignmentRecord],
) -> io::Result<(AdaptiveCodec, AdaptiveCodec)> {
    let mut first_val_byte_frequencies: HashMap<i64, u64> = HashMap::new();
    let mut second_val_byte_frequencies: HashMap<i64, u64> = HashMap::new();

    // Collect byte frequencies by varint encoding all raw values
    for record in records {
        let tracepoints = match &record.tracepoints {
            TracepointData::Standard(tps) | TracepointData::Fastga(tps) => tps,
            _ => continue,
        };

        if tracepoints.is_empty() {
            continue;
        }

        // Encode all values as raw (no delta)
        for &(first_val, second_val) in tracepoints {
            // First value as raw
            let mut buf = Vec::new();
            write_varint(&mut buf, first_val).unwrap();
            for byte in buf {
                *first_val_byte_frequencies.entry(byte as i64).or_insert(0) += 1;
            }

            // Second value as raw
            let mut buf = Vec::new();
            write_varint(&mut buf, second_val).unwrap();
            for byte in buf {
                *second_val_byte_frequencies.entry(byte as i64).or_insert(0) += 1;
            }
        }
    }

    if first_val_byte_frequencies.is_empty() || second_val_byte_frequencies.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "No tracepoint data found for training byte codecs",
        ));
    }

    // Debug: log byte frequency distribution statistics
    debug!("First value byte statistics (raw, after varint):");
    debug!("  Unique byte symbols: {}", first_val_byte_frequencies.len());
    let first_total: u64 = first_val_byte_frequencies.values().sum();
    let mut first_sorted: Vec<_> = first_val_byte_frequencies.iter().collect();
    first_sorted.sort_by(|a, b| b.1.cmp(a.1));
    if let Some((val, freq)) = first_sorted.first() {
        debug!("  Most frequent: 0x{:02X} ({}x, {:.2}%)", val, freq, (**freq as f64 / first_total as f64) * 100.0);
    }
    if first_sorted.len() >= 3 {
        let top3_freq: u64 = first_sorted.iter().take(3).map(|(_, f)| **f).sum();
        debug!("  Top 3 bytes: {:.2}% of data", (top3_freq as f64 / first_total as f64) * 100.0);
    }

    debug!("Second value byte statistics (raw, after varint):");
    debug!("  Unique byte symbols: {}", second_val_byte_frequencies.len());
    let second_total: u64 = second_val_byte_frequencies.values().sum();
    let mut second_sorted: Vec<_> = second_val_byte_frequencies.iter().collect();
    second_sorted.sort_by(|a, b| b.1.cmp(a.1));
    if let Some((val, freq)) = second_sorted.first() {
        debug!("  Most frequent: 0x{:02X} ({}x, {:.2}%)", val, freq, (**freq as f64 / second_total as f64) * 100.0);
    }
    if second_sorted.len() >= 3 {
        let top3_freq: u64 = second_sorted.iter().take(3).map(|(_, f)| **f).sum();
        debug!("  Top 3 bytes: {:.2}% of data", (top3_freq as f64 / second_total as f64) * 100.0);
    }

    let first_val_codec = AdaptiveCodec::build(&first_val_byte_frequencies)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    let second_val_codec = AdaptiveCodec::build(&second_val_byte_frequencies)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    Ok((first_val_codec, second_val_codec))
}

/// Analyze data distribution and select optimal compression strategy
///
/// Tests all 6 strategies by actually compressing the data and measuring size
/// Result from Smart strategy selection
struct SmartStrategyResult {
    strategy: CompressionStrategy,
    codecs: Option<(AdaptiveCodec, AdaptiveCodec)>,
}

/// Returns the strategy that produces the smallest compressed output
/// For Huffman-based strategies, also returns the trained codecs to avoid re-training
fn analyze_and_select_strategy(
    records: &[AlignmentRecord],
) -> io::Result<SmartStrategyResult> {
    use std::io::Cursor;

    info!("Smart: Testing all 6 compression strategies...");

    // Test each strategy by encoding to memory
    let strategies = vec![
        CompressionStrategy::Varint,
        CompressionStrategy::VarintRaw,
        CompressionStrategy::Huffman,
        CompressionStrategy::HuffmanRaw,
        CompressionStrategy::VarintHuffman,
        CompressionStrategy::VarintHuffmanRaw,
    ];

    let mut results: Vec<(CompressionStrategy, usize, Option<(AdaptiveCodec, AdaptiveCodec)>)> = Vec::new();

    for strategy in strategies {
        // Create a temporary in-memory buffer
        let mut buffer = Cursor::new(Vec::new());

        // Encode tracepoints only (the main compressible data)
        // For Huffman strategies, also capture the trained codecs
        let result: io::Result<(usize, Option<(AdaptiveCodec, AdaptiveCodec)>)> = match strategy {
            CompressionStrategy::Varint | CompressionStrategy::VarintRaw => {
                // Encode with varint (matching real writer: separate first/second compression)
                for record in records {
                    let tracepoints = match &record.tracepoints {
                        TracepointData::Standard(tps) | TracepointData::Fastga(tps) => tps,
                        _ => continue,
                    };

                    if tracepoints.is_empty() {
                        continue;
                    }

                    let first_vals: Vec<u64> = tracepoints.iter().map(|&(v, _)| v).collect();

                    // Encode first values
                    let mut first_val_buf = Vec::new();
                    if strategy == CompressionStrategy::Varint {
                        // Varint: FastGA-aware delta encoding with zigzag
                        let use_delta = !matches!(record.tp_type, TracepointType::Fastga);
                        let first_vals_encoded: Vec<i64> = if use_delta {
                            delta_encode(&first_vals)
                        } else {
                            first_vals.iter().map(|&v| v as i64).collect()
                        };
                        for val in first_vals_encoded {
                            let zigzag = ((val << 1) ^ (val >> 63)) as u64;
                            write_varint(&mut first_val_buf, zigzag)?;
                        }
                    } else {
                        // VarintRaw: raw values, no zigzag
                        for &val in &first_vals {
                            write_varint(&mut first_val_buf, val)?;
                        }
                    }
                    let first_compressed = zstd::encode_all(&first_val_buf[..], 3)?;
                    write_varint(&mut buffer, first_compressed.len() as u64)?;
                    buffer.write_all(&first_compressed)?;

                    // Encode second values
                    let mut second_val_buf = Vec::new();
                    for &(_, val) in tracepoints {
                        write_varint(&mut second_val_buf, val)?;
                    }
                    let second_compressed = zstd::encode_all(&second_val_buf[..], 3)?;
                    write_varint(&mut buffer, second_compressed.len() as u64)?;
                    buffer.write_all(&second_compressed)?;
                }
                Ok((buffer.into_inner().len(), None))
            }
            CompressionStrategy::Huffman | CompressionStrategy::HuffmanRaw => {
                // Train and encode with integer-level Huffman
                let (first_codec, second_codec) = if strategy == CompressionStrategy::Huffman {
                    train_codecs_from_records(records)?
                } else {
                    train_codecs_from_records_raw(records)?
                };

                // Include codec serialization overhead in size calculation
                first_codec.write(&mut buffer)?;
                second_codec.write(&mut buffer)?;

                for record in records {
                    let tracepoints = match &record.tracepoints {
                        TracepointData::Standard(tps) | TracepointData::Fastga(tps) => tps,
                        _ => continue,
                    };

                    let use_delta = strategy == CompressionStrategy::Huffman && !matches!(record.tp_type, TracepointType::Fastga);

                    let first_vals: Vec<u64> = tracepoints.iter().map(|&(v, _)| v).collect();
                    let first_vals_encoded: Vec<i64> = if use_delta {
                        delta_encode(&first_vals)
                    } else {
                        first_vals.iter().map(|&v| v as i64).collect()
                    };

                    let mut first_writer = BitWriter::new();
                    for val in first_vals_encoded {
                        let code = first_codec.encode(val).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                        first_writer.write_bits(&code.bits);
                    }
                    let first_bytes = first_writer.finish();
                    let first_compressed = zstd::encode_all(&first_bytes[..], 3)?;
                    write_varint(&mut buffer, first_compressed.len() as u64)?;
                    buffer.write_all(&first_compressed)?;

                    let mut second_writer = BitWriter::new();
                    for &(_, val) in tracepoints {
                        let code = second_codec.encode(val as i64).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                        second_writer.write_bits(&code.bits);
                    }
                    let second_bytes = second_writer.finish();
                    let second_compressed = zstd::encode_all(&second_bytes[..], 3)?;
                    write_varint(&mut buffer, second_compressed.len() as u64)?;
                    buffer.write_all(&second_compressed)?;
                }
                Ok((buffer.into_inner().len(), Some((first_codec, second_codec))))
            }
            CompressionStrategy::VarintHuffman | CompressionStrategy::VarintHuffmanRaw => {
                // Train and encode with byte-level Huffman
                let (first_codec, second_codec) = if strategy == CompressionStrategy::VarintHuffman {
                    train_codecs_from_bytes(records)?
                } else {
                    train_codecs_from_bytes_raw(records)?
                };

                // Include codec serialization overhead in size calculation
                first_codec.write(&mut buffer)?;
                second_codec.write(&mut buffer)?;

                for record in records {
                    let tracepoints = match &record.tracepoints {
                        TracepointData::Standard(tps) | TracepointData::Fastga(tps) => tps,
                        _ => continue,
                    };

                    let use_delta = strategy == CompressionStrategy::VarintHuffman && !matches!(record.tp_type, TracepointType::Fastga);

                    let first_vals: Vec<u64> = tracepoints.iter().map(|&(v, _)| v).collect();

                    let mut first_bytes_buf = Vec::new();
                    if strategy == CompressionStrategy::VarintHuffman {
                        // VarintHuffman: use zigzag encoding (for both delta and raw)
                        let first_vals_encoded: Vec<i64> = if use_delta {
                            delta_encode(&first_vals)
                        } else {
                            first_vals.iter().map(|&v| v as i64).collect()
                        };
                        for val in first_vals_encoded {
                            let zigzag = ((val << 1) ^ (val >> 63)) as u64;
                            write_varint(&mut first_bytes_buf, zigzag)?;
                        }
                    } else {
                        // VarintHuffmanRaw: no zigzag, plain varint on raw u64
                        for &val in &first_vals {
                            write_varint(&mut first_bytes_buf, val)?;
                        }
                    }

                    let mut first_writer = BitWriter::new();
                    for byte in first_bytes_buf {
                        let code = first_codec.encode(byte as i64).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                        first_writer.write_bits(&code.bits);
                    }
                    let first_bytes = first_writer.finish();
                    let first_compressed = zstd::encode_all(&first_bytes[..], 3)?;
                    write_varint(&mut buffer, first_compressed.len() as u64)?;
                    buffer.write_all(&first_compressed)?;

                    let mut second_bytes_buf = Vec::new();
                    for &(_, val) in tracepoints {
                        write_varint(&mut second_bytes_buf, val)?;
                    }

                    let mut second_writer = BitWriter::new();
                    for byte in second_bytes_buf {
                        let code = second_codec.encode(byte as i64).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                        second_writer.write_bits(&code.bits);
                    }
                    let second_bytes = second_writer.finish();
                    let second_compressed = zstd::encode_all(&second_bytes[..], 3)?;
                    write_varint(&mut buffer, second_compressed.len() as u64)?;
                    buffer.write_all(&second_compressed)?;
                }
                Ok((buffer.into_inner().len(), Some((first_codec, second_codec))))
            }
            CompressionStrategy::Smart => {
                unreachable!("Smart should not test itself")
            }
        };

        match result {
            Ok((size, codecs)) => {
                info!("Smart: {} → {} bytes", strategy, size);
                results.push((strategy, size, codecs));
            }
            Err(e) => {
                debug!("Smart: {} failed: {}", strategy, e);
            }
        }
    }

    // Find the strategy with minimum size
    if let Some((best_strategy, best_size, best_codecs)) = results.iter().min_by_key(|(_, size, _)| size) {
        info!("Smart: Selected {} ({} bytes, best compression)", best_strategy, best_size);
        Ok(SmartStrategyResult {
            strategy: *best_strategy,
            codecs: best_codecs.clone(),
        })
    } else {
        // Fallback to Varint if all strategies failed
        info!("Smart: All strategies failed, defaulting to Varint");
        Ok(SmartStrategyResult {
            strategy: CompressionStrategy::Varint,
            codecs: None,
        })
    }
}

/// Write records to binary PAF format with Huffman encoding
fn write_binary_adaptive(
    output_path: &str,
    records: &[AlignmentRecord],
    string_table: &StringTable,
    first_val_codec: &AdaptiveCodec,
    second_val_codec: &AdaptiveCodec,
) -> io::Result<()> {
    let output = File::create(output_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to create output file '{}': {}", output_path, e),
        )
    })?;
    let mut writer = BufWriter::new(output);

    let header = BinaryPafHeader {
        version: 1,
        flags: CompressionStrategy::Huffman.to_code(),
        num_records: records.len() as u64,
        num_strings: string_table.len() as u64,
    };

    header.write(&mut writer)?;

    // Write codecs
    first_val_codec.write(&mut writer)?;
    second_val_codec.write(&mut writer)?;

    // Write records
    for record in records {
        write_record_adaptive(&mut writer, record, first_val_codec, second_val_codec)?;
    }

    // Write string table
    string_table.write(&mut writer)?;

    writer.flush()?;
    Ok(())
}

/// Write single record with Huffman encoding
fn write_record_adaptive<W: Write>(
    writer: &mut W,
    record: &AlignmentRecord,
    first_val_codec: &AdaptiveCodec,
    second_val_codec: &AdaptiveCodec,
) -> io::Result<()> {
    // Write core PAF fields
    write_varint(writer, record.query_name_id)?;
    write_varint(writer, record.query_start)?;
    write_varint(writer, record.query_end)?;
    writer.write_all(&[record.strand as u8])?;
    write_varint(writer, record.target_name_id)?;
    write_varint(writer, record.target_start)?;
    write_varint(writer, record.target_end)?;
    write_varint(writer, record.residue_matches)?;
    write_varint(writer, record.alignment_block_len)?;
    writer.write_all(&[record.mapping_quality])?;
    writer.write_all(&[record.tp_type.to_u8()])?;
    writer.write_all(&[complexity_metric_to_u8(&record.complexity_metric)])?;
    write_varint(writer, record.max_complexity)?;

    // Get tracepoints
    let empty_vec = Vec::new();
    let tracepoints = match &record.tracepoints {
        TracepointData::Standard(tps) | TracepointData::Fastga(tps) => tps,
        _ => &empty_vec,
    };

    write_varint(writer, tracepoints.len() as u64)?;

    if !tracepoints.is_empty() {
        // FastGA-aware delta encoding:
        // - FastGA: num_differences are naturally small, use raw values
        // - Standard: query_bases are incremental, use delta encoding
        let use_delta = !matches!(record.tp_type, TracepointType::Fastga);

        // Encode first values with Huffman
        let mut first_val_writer = BitWriter::new();
        let first_vals: Vec<u64> = tracepoints.iter().map(|&(v, _)| v).collect();
        let first_vals_encoded: Vec<i64> = if use_delta {
            delta_encode(&first_vals)
        } else {
            first_vals.iter().map(|&v| v as i64).collect()
        };

        for val in first_vals_encoded {
            let code = first_val_codec
                .encode(val)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            first_val_writer.write_bits(&code.bits);
        }
        let first_val_bytes = first_val_writer.finish();
        let first_val_compressed = zstd::encode_all(&first_val_bytes[..], 3)?;
        write_varint(writer, first_val_compressed.len() as u64)?;
        writer.write_all(&first_val_compressed)?;

        // Encode second values with Huffman
        let mut second_val_writer = BitWriter::new();
        for &(_, val) in tracepoints {
            let code = second_val_codec
                .encode(val as i64)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            second_val_writer.write_bits(&code.bits);
        }
        let second_val_bytes = second_val_writer.finish();
        let second_val_compressed = zstd::encode_all(&second_val_bytes[..], 3)?;
        write_varint(writer, second_val_compressed.len() as u64)?;
        writer.write_all(&second_val_compressed)?;
    }

    // Write tags
    write_varint(writer, record.tags.len() as u64)?;
    for tag in &record.tags {
        tag.write(writer)?;
    }

    Ok(())
}

/// Write single record with Huffman encoding (NO delta encoding)
fn write_record_adaptive_raw<W: Write>(
    writer: &mut W,
    record: &AlignmentRecord,
    first_val_codec: &AdaptiveCodec,
    second_val_codec: &AdaptiveCodec,
) -> io::Result<()> {
    // Write core PAF fields
    write_varint(writer, record.query_name_id)?;
    write_varint(writer, record.query_start)?;
    write_varint(writer, record.query_end)?;
    writer.write_all(&[record.strand as u8])?;
    write_varint(writer, record.target_name_id)?;
    write_varint(writer, record.target_start)?;
    write_varint(writer, record.target_end)?;
    write_varint(writer, record.residue_matches)?;
    write_varint(writer, record.alignment_block_len)?;
    writer.write_all(&[record.mapping_quality])?;
    writer.write_all(&[record.tp_type.to_u8()])?;
    writer.write_all(&[complexity_metric_to_u8(&record.complexity_metric)])?;
    write_varint(writer, record.max_complexity)?;

    // Get tracepoints
    let empty_vec = Vec::new();
    let tracepoints = match &record.tracepoints {
        TracepointData::Standard(tps) | TracepointData::Fastga(tps) => tps,
        _ => &empty_vec,
    };

    write_varint(writer, tracepoints.len() as u64)?;

    if !tracepoints.is_empty() {
        // No delta encoding - use raw values for all types
        let mut first_val_writer = BitWriter::new();
        for &(val, _) in tracepoints {
            let code = first_val_codec
                .encode(val as i64)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            first_val_writer.write_bits(&code.bits);
        }
        let first_val_bytes = first_val_writer.finish();
        let first_val_compressed = zstd::encode_all(&first_val_bytes[..], 3)?;
        write_varint(writer, first_val_compressed.len() as u64)?;
        writer.write_all(&first_val_compressed)?;

        // Encode second values with Huffman
        let mut second_val_writer = BitWriter::new();
        for &(_, val) in tracepoints {
            let code = second_val_codec
                .encode(val as i64)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            second_val_writer.write_bits(&code.bits);
        }
        let second_val_bytes = second_val_writer.finish();
        let second_val_compressed = zstd::encode_all(&second_val_bytes[..], 3)?;
        write_varint(writer, second_val_compressed.len() as u64)?;
        writer.write_all(&second_val_compressed)?;
    }

    // Write tags
    write_varint(writer, record.tags.len() as u64)?;
    for tag in &record.tags {
        tag.write(writer)?;
    }

    Ok(())
}

/// Write binary with Huffman encoding (NO delta encoding)
fn write_binary_adaptive_raw(
    output_path: &str,
    records: &[AlignmentRecord],
    string_table: &StringTable,
    first_val_codec: &AdaptiveCodec,
    second_val_codec: &AdaptiveCodec,
) -> io::Result<()> {
    let output = File::create(output_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to create output file '{}': {}", output_path, e),
        )
    })?;
    let mut writer = BufWriter::new(output);

    let header = BinaryPafHeader {
        version: 1,
        flags: CompressionStrategy::HuffmanRaw.to_code(),
        num_records: records.len() as u64,
        num_strings: string_table.len() as u64,
    };

    header.write(&mut writer)?;

    // Write codecs
    first_val_codec.write(&mut writer)?;
    second_val_codec.write(&mut writer)?;

    // Write records
    for record in records {
        write_record_adaptive_raw(&mut writer, record, first_val_codec, second_val_codec)?;
    }

    // Write string table
    string_table.write(&mut writer)?;

    writer.flush()?;
    Ok(())
}

/// Write records with hybrid varint+Huffman compression
///
/// Pipeline: Delta/Raw → Varint → Byte-level Huffman → Zstd
/// Format version 5 = varint-huffman hybrid
fn write_binary_varint_huffman(
    output_path: &str,
    records: &[AlignmentRecord],
    string_table: &StringTable,
    first_val_codec: &AdaptiveCodec,
    second_val_codec: &AdaptiveCodec,
) -> io::Result<()> {
    let output = File::create(output_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to create output file '{}': {}", output_path, e),
        )
    })?;
    let mut writer = BufWriter::new(output);

    let header = BinaryPafHeader {
        version: 1,
        flags: CompressionStrategy::VarintHuffman.to_code(),
        num_records: records.len() as u64,
        num_strings: string_table.len() as u64,
    };

    header.write(&mut writer)?;

    // Write codecs
    first_val_codec.write(&mut writer)?;
    second_val_codec.write(&mut writer)?;

    // Write records
    for record in records {
        write_record_varint_huffman(&mut writer, record, first_val_codec, second_val_codec)?;
    }

    // Write string table
    string_table.write(&mut writer)?;

    writer.flush()?;
    Ok(())
}

/// Write single record with varint+Huffman hybrid encoding
fn write_record_varint_huffman<W: Write>(
    writer: &mut W,
    record: &AlignmentRecord,
    first_val_codec: &AdaptiveCodec,
    second_val_codec: &AdaptiveCodec,
) -> io::Result<()> {
    // Write core PAF fields
    write_varint(writer, record.query_name_id)?;
    write_varint(writer, record.query_start)?;
    write_varint(writer, record.query_end)?;
    writer.write_all(&[record.strand as u8])?;
    write_varint(writer, record.target_name_id)?;
    write_varint(writer, record.target_start)?;
    write_varint(writer, record.target_end)?;
    write_varint(writer, record.residue_matches)?;
    write_varint(writer, record.alignment_block_len)?;
    writer.write_all(&[record.mapping_quality])?;
    writer.write_all(&[record.tp_type.to_u8()])?;
    writer.write_all(&[complexity_metric_to_u8(&record.complexity_metric)])?;
    write_varint(writer, record.max_complexity)?;

    // Get tracepoints
    let empty_vec = Vec::new();
    let tracepoints = match &record.tracepoints {
        TracepointData::Standard(tps) | TracepointData::Fastga(tps) => tps,
        _ => &empty_vec,
    };

    write_varint(writer, tracepoints.len() as u64)?;

    if !tracepoints.is_empty() {
        // FastGA-aware delta encoding
        let use_delta = !matches!(record.tp_type, TracepointType::Fastga);

        // Encode first values: delta/raw → varint → bytes
        let first_vals: Vec<u64> = tracepoints.iter().map(|&(v, _)| v).collect();
        let first_vals_encoded: Vec<i64> = if use_delta {
            delta_encode(&first_vals)
        } else {
            first_vals.iter().map(|&v| v as i64).collect()
        };

        // Varint encode to bytes
        let mut first_val_bytes_buf = Vec::new();
        for val in first_vals_encoded {
            let zigzag = ((val << 1) ^ (val >> 63)) as u64;
            write_varint(&mut first_val_bytes_buf, zigzag)?;
        }

        // Apply byte-level Huffman on varint output
        let mut first_val_writer = BitWriter::new();
        for byte in first_val_bytes_buf {
            let code = first_val_codec
                .encode(byte as i64)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            first_val_writer.write_bits(&code.bits);
        }
        let first_val_bytes = first_val_writer.finish();
        let first_val_compressed = zstd::encode_all(&first_val_bytes[..], 3)?;
        write_varint(writer, first_val_compressed.len() as u64)?;
        writer.write_all(&first_val_compressed)?;

        // Encode second values: raw → varint → bytes
        let mut second_val_bytes_buf = Vec::new();
        for &(_, val) in tracepoints {
            write_varint(&mut second_val_bytes_buf, val)?;
        }

        // Apply byte-level Huffman on varint output
        let mut second_val_writer = BitWriter::new();
        for byte in second_val_bytes_buf {
            let code = second_val_codec
                .encode(byte as i64)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            second_val_writer.write_bits(&code.bits);
        }
        let second_val_bytes = second_val_writer.finish();
        let second_val_compressed = zstd::encode_all(&second_val_bytes[..], 3)?;
        write_varint(writer, second_val_compressed.len() as u64)?;
        writer.write_all(&second_val_compressed)?;
    }

    // Write tags
    write_varint(writer, record.tags.len() as u64)?;
    for tag in &record.tags {
        tag.write(writer)?;
    }

    Ok(())
}

/// Write records with hybrid varint+Huffman compression (NO delta encoding)
///
/// Pipeline: Raw → Varint → Byte-level Huffman → Zstd
/// Format version 6 = varint-huffman-raw hybrid
fn write_binary_varint_huffman_raw(
    output_path: &str,
    records: &[AlignmentRecord],
    string_table: &StringTable,
    first_val_codec: &AdaptiveCodec,
    second_val_codec: &AdaptiveCodec,
) -> io::Result<()> {
    let output = File::create(output_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to create output file '{}': {}", output_path, e),
        )
    })?;
    let mut writer = BufWriter::new(output);

    let header = BinaryPafHeader {
        version: 1,
        flags: CompressionStrategy::VarintHuffmanRaw.to_code(),
        num_records: records.len() as u64,
        num_strings: string_table.len() as u64,
    };

    header.write(&mut writer)?;

    // Write codecs
    first_val_codec.write(&mut writer)?;
    second_val_codec.write(&mut writer)?;

    // Write records
    for record in records {
        write_record_varint_huffman_raw(&mut writer, record, first_val_codec, second_val_codec)?;
    }

    // Write string table
    string_table.write(&mut writer)?;

    writer.flush()?;
    Ok(())
}

/// Write single record with varint+Huffman hybrid encoding (NO delta)
fn write_record_varint_huffman_raw<W: Write>(
    writer: &mut W,
    record: &AlignmentRecord,
    first_val_codec: &AdaptiveCodec,
    second_val_codec: &AdaptiveCodec,
) -> io::Result<()> {
    // Write core PAF fields
    write_varint(writer, record.query_name_id)?;
    write_varint(writer, record.query_start)?;
    write_varint(writer, record.query_end)?;
    writer.write_all(&[record.strand as u8])?;
    write_varint(writer, record.target_name_id)?;
    write_varint(writer, record.target_start)?;
    write_varint(writer, record.target_end)?;
    write_varint(writer, record.residue_matches)?;
    write_varint(writer, record.alignment_block_len)?;
    writer.write_all(&[record.mapping_quality])?;
    writer.write_all(&[record.tp_type.to_u8()])?;
    writer.write_all(&[complexity_metric_to_u8(&record.complexity_metric)])?;
    write_varint(writer, record.max_complexity)?;

    // Get tracepoints
    let empty_vec = Vec::new();
    let tracepoints = match &record.tracepoints {
        TracepointData::Standard(tps) | TracepointData::Fastga(tps) => tps,
        _ => &empty_vec,
    };

    write_varint(writer, tracepoints.len() as u64)?;

    if !tracepoints.is_empty() {
        // Encode first values: raw → varint → bytes (NO delta)
        let mut first_val_bytes_buf = Vec::new();
        for &(val, _) in tracepoints {
            write_varint(&mut first_val_bytes_buf, val)?;
        }

        // Apply byte-level Huffman on varint output
        let mut first_val_writer = BitWriter::new();
        for byte in first_val_bytes_buf {
            let code = first_val_codec
                .encode(byte as i64)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            first_val_writer.write_bits(&code.bits);
        }
        let first_val_bytes = first_val_writer.finish();
        let first_val_compressed = zstd::encode_all(&first_val_bytes[..], 3)?;
        write_varint(writer, first_val_compressed.len() as u64)?;
        writer.write_all(&first_val_compressed)?;

        // Encode second values: raw → varint → bytes
        let mut second_val_bytes_buf = Vec::new();
        for &(_, val) in tracepoints {
            write_varint(&mut second_val_bytes_buf, val)?;
        }

        // Apply byte-level Huffman on varint output
        let mut second_val_writer = BitWriter::new();
        for byte in second_val_bytes_buf {
            let code = second_val_codec
                .encode(byte as i64)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            second_val_writer.write_bits(&code.bits);
        }
        let second_val_bytes = second_val_writer.finish();
        let second_val_compressed = zstd::encode_all(&second_val_bytes[..], 3)?;
        write_varint(writer, second_val_compressed.len() as u64)?;
        writer.write_all(&second_val_compressed)?;
    }

    // Write tags
    write_varint(writer, record.tags.len() as u64)?;
    for tag in &record.tags {
        tag.write(writer)?;
    }

    Ok(())
}

/// Write records to binary PAF format
///
/// Uses FastGA-aware delta encoding, varint, and zstd compression for efficient storage.
fn write_binary(
    output_path: &str,
    records: &[AlignmentRecord],
    string_table: &StringTable,
) -> io::Result<()> {
    let output = File::create(output_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to create output file '{}': {}", output_path, e),
        )
    })?;
    let mut writer = BufWriter::new(output);

    let header = BinaryPafHeader {
        version: 1,
        flags: CompressionStrategy::Varint.to_code(),
        num_records: records.len() as u64,
        num_strings: string_table.len() as u64,
    };

    header.write(&mut writer)?;
    for record in records {
        record.write(&mut writer)?;
    }
    string_table.write(&mut writer)?;
    Ok(())
}

/// Write records to binary PAF format (NO delta encoding)
///
/// Uses raw values, varint, and zstd compression.
fn write_binary_raw(
    output_path: &str,
    records: &[AlignmentRecord],
    string_table: &StringTable,
) -> io::Result<()> {
    let output = File::create(output_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("Failed to create output file '{}': {}", output_path, e),
        )
    })?;
    let mut writer = BufWriter::new(output);

    let header = BinaryPafHeader {
        version: 1,
        flags: CompressionStrategy::VarintRaw.to_code(),
        num_records: records.len() as u64,
        num_strings: string_table.len() as u64,
    };

    header.write(&mut writer)?;
    for record in records {
        record.write_raw(&mut writer)?;
    }
    string_table.write(&mut writer)?;
    Ok(())
}

fn write_paf_output(
    output_path: &str,
    records: &[AlignmentRecord],
    string_table: &StringTable,
) -> io::Result<()> {
    let output: Box<dyn Write> = if output_path == "-" {
        Box::new(io::stdout())
    } else {
        Box::new(File::create(output_path).map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("Failed to create output file '{}': {}", output_path, e),
            )
        })?)
    };
    let mut writer = BufWriter::new(output);

    for record in records {
        write_paf_line(&mut writer, record, string_table)?;
    }
    writer.flush()?;
    Ok(())
}

fn write_paf_line<W: Write>(
    writer: &mut W,
    record: &AlignmentRecord,
    string_table: &StringTable,
) -> io::Result<()> {
    let query_name = string_table.get(record.query_name_id).unwrap();
    let target_name = string_table.get(record.target_name_id).unwrap();
    let query_len = string_table.get_length(record.query_name_id).unwrap();
    let target_len = string_table.get_length(record.target_name_id).unwrap();

    write!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        query_name,
        query_len,
        record.query_start,
        record.query_end,
        record.strand,
        target_name,
        target_len,
        record.target_start,
        record.target_end,
        record.residue_matches,
        record.alignment_block_len,
        record.mapping_quality
    )?;

    for tag in &record.tags {
        write!(writer, "\t{}", format_tag(tag))?;
    }
    write!(writer, "\ttp:Z:{}", format_tracepoints(&record.tracepoints))?;
    writeln!(writer)?;
    Ok(())
}

fn format_tracepoints(tps: &TracepointData) -> String {
    match tps {
        TracepointData::Standard(items) | TracepointData::Fastga(items) => items
            .iter()
            .map(|(a, b)| format!("{},{}", a, b))
            .collect::<Vec<_>>()
            .join(";"),
        TracepointData::Variable(items) => items
            .iter()
            .map(|(a, b_opt)| match b_opt {
                Some(b) => format!("{},{}", a, b),
                None => a.to_string(),
            })
            .collect::<Vec<_>>()
            .join(";"),
        TracepointData::Mixed(items) => items
            .iter()
            .map(|item| match item {
                MixedTracepointItem::Tracepoint(a, b) => format!("{},{}", a, b),
                MixedTracepointItem::CigarOp(len, op) => format!("{}{}", len, op),
            })
            .collect::<Vec<_>>()
            .join(";"),
    }
}

fn parse_tag(field: &str) -> Option<Tag> {
    let parts: Vec<&str> = field.splitn(3, ':').collect();
    if parts.len() != 3 {
        return None;
    }
    let key = [parts[0].as_bytes()[0], parts[0].as_bytes()[1]];
    let tag_type = parts[1].as_bytes()[0];
    let value = match tag_type {
        b'i' => parts[2].parse::<i64>().ok().map(TagValue::Int)?,
        b'f' => parts[2].parse::<f32>().ok().map(TagValue::Float)?,
        b'Z' => TagValue::String(parts[2].to_string()),
        _ => return None,
    };
    Some(Tag {
        key,
        tag_type,
        value,
    })
}

fn format_tag(tag: &Tag) -> String {
    let key = String::from_utf8_lossy(&tag.key);
    match &tag.value {
        TagValue::Int(v) => format!("{}:i:{}", key, v),
        TagValue::Float(v) => format!("{}:f:{}", key, v),
        TagValue::String(s) => format!("{}:Z:{}", key, s),
    }
}

fn parse_usize(s: &str, field: &str) -> io::Result<u64> {
    s.parse().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid {}: {}", field, s),
        )
    })
}

fn parse_u8(s: &str, field: &str) -> io::Result<u8> {
    s.parse().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid {}: {}", field, s),
        )
    })
}

/// Convert ComplexityMetric to u8 for serialization
#[inline]
fn complexity_metric_to_u8(metric: &ComplexityMetric) -> u8 {
    match metric {
        ComplexityMetric::EditDistance => 0,
        ComplexityMetric::DiagonalDistance => 1,
    }
}

/// Convert u8 to ComplexityMetric for deserialization
#[inline]
fn complexity_metric_from_u8(byte: u8) -> io::Result<ComplexityMetric> {
    match byte {
        0 => Ok(ComplexityMetric::EditDistance),
        1 => Ok(ComplexityMetric::DiagonalDistance),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Invalid complexity metric",
        )),
    }
}

fn parse_paf_with_cigar(
    line: &str,
    string_table: &mut StringTable,
    tp_type: &TracepointType,
    max_complexity: usize,
    complexity_metric: &ComplexityMetric,
) -> io::Result<AlignmentRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "PAF line has fewer than 12 fields",
        ));
    }

    let query_len = parse_usize(fields[1], "query_len")?;
    let query_name_id = string_table.intern(fields[0], query_len);
    let query_start = parse_usize(fields[2], "query_start")?;
    let query_end = parse_usize(fields[3], "query_end")?;
    let strand = fields[4].chars().next().unwrap_or('+');
    let target_len = parse_usize(fields[6], "target_len")?;
    let target_name_id = string_table.intern(fields[5], target_len);
    let target_start = parse_usize(fields[7], "target_start")?;
    let target_end = parse_usize(fields[8], "target_end")?;
    let residue_matches = parse_usize(fields[9], "residue_matches")?;
    let alignment_block_len = parse_usize(fields[10], "alignment_block_len")?;
    let mapping_quality = parse_u8(fields[11], "mapping_quality")?;

    let mut cigar = None;
    let mut tags = Vec::new();

    for field in &fields[12..] {
        if let Some(stripped) = field.strip_prefix("cg:Z:") {
            cigar = Some(stripped);
        } else if !field.starts_with("tp:Z:") {
            if let Some(tag) = parse_tag(field) {
                tags.push(tag);
            }
        }
    }

    let cigar =
        cigar.ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing cg:Z: tag"))?;

    let tracepoints = match tp_type {
        TracepointType::Standard => {
            let tps = cigar_to_tracepoints(cigar, max_complexity, *complexity_metric);
            TracepointData::Standard(tps.into_iter().map(|(a, b)| (a as u64, b as u64)).collect())
        }
        TracepointType::Mixed => {
            let tps = cigar_to_mixed_tracepoints(cigar, max_complexity, *complexity_metric);
            TracepointData::Mixed(tps.iter().map(MixedTracepointItem::from).collect())
        }
        TracepointType::Variable => {
            let tps = cigar_to_variable_tracepoints(cigar, max_complexity, *complexity_metric);
            TracepointData::Variable(
                tps.into_iter()
                    .map(|(a, b_opt)| (a as u64, b_opt.map(|b| b as u64)))
                    .collect(),
            )
        }
        TracepointType::Fastga => {
            let complement = strand == '-';
            let segments = cigar_to_tracepoints_fastga(
                cigar,
                max_complexity,
                query_start as usize,
                query_end as usize,
                query_len as usize,
                target_start as usize,
                target_end as usize,
                target_len as usize,
                complement,
            );
            if let Some((tps, _coords)) = segments.first() {
                TracepointData::Fastga(tps.iter().map(|(a, b)| (*a as u64, *b as u64)).collect())
            } else {
                TracepointData::Fastga(Vec::new())
            }
        }
    };

    let tp_type_enum = *tp_type;

    Ok(AlignmentRecord {
        query_name_id,
        query_start,
        query_end,
        strand,
        target_name_id,
        target_start,
        target_end,
        residue_matches,
        alignment_block_len,
        mapping_quality,
        tp_type: tp_type_enum,
        complexity_metric: *complexity_metric,
        max_complexity: max_complexity as u64,
        tracepoints,
        tags,
    })
}

fn parse_paf_with_tracepoints(
    line: &str,
    string_table: &mut StringTable,
) -> io::Result<AlignmentRecord> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "PAF line has fewer than 12 fields",
        ));
    }

    let query_len = parse_usize(fields[1], "query_len")?;
    let query_name_id = string_table.intern(fields[0], query_len);
    let query_start = parse_usize(fields[2], "query_start")?;
    let query_end = parse_usize(fields[3], "query_end")?;
    let strand = fields[4].chars().next().unwrap_or('+');
    let target_len = parse_usize(fields[6], "target_len")?;
    let target_name_id = string_table.intern(fields[5], target_len);
    let target_start = parse_usize(fields[7], "target_start")?;
    let target_end = parse_usize(fields[8], "target_end")?;
    let residue_matches = parse_usize(fields[9], "residue_matches")?;
    let alignment_block_len = parse_usize(fields[10], "alignment_block_len")?;
    let mapping_quality = parse_u8(fields[11], "mapping_quality")?;

    let mut tp_str = None;
    let mut tags = Vec::new();

    for field in &fields[12..] {
        if let Some(stripped) = field.strip_prefix("tp:Z:") {
            tp_str = Some(stripped);
        } else if let Some(tag) = parse_tag(field) {
            tags.push(tag);
        }
    }

    let tp_str =
        tp_str.ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Missing tp:Z: tag"))?;
    let (tracepoints, tp_type) = parse_tracepoints_auto(tp_str)?;

    Ok(AlignmentRecord {
        query_name_id,
        query_start,
        query_end,
        strand,
        target_name_id,
        target_start,
        target_end,
        residue_matches,
        alignment_block_len,
        mapping_quality,
        tp_type,
        complexity_metric: ComplexityMetric::EditDistance,
        max_complexity: 100,
        tracepoints,
        tags,
    })
}

fn parse_tracepoints_auto(tp_str: &str) -> io::Result<(TracepointData, TracepointType)> {
    if tp_str.is_empty() {
        return Ok((
            TracepointData::Standard(Vec::new()),
            TracepointType::Standard,
        ));
    }

    let items: Vec<&str> = tp_str.split(';').collect();
    let has_cigar = items.iter().any(|s| {
        s.chars()
            .last()
            .map(|c| matches!(c, 'M' | 'D' | 'I' | 'X' | '='))
            .unwrap_or(false)
    });

    if has_cigar {
        let mut mixed = Vec::new();
        for item in items {
            if item
                .chars()
                .last()
                .map(|c| matches!(c, 'M' | 'D' | 'I' | 'X' | '='))
                .unwrap_or(false)
            {
                let op = item.chars().last().unwrap();
                let len_str = &item[..item.len() - 1];
                let len: u64 = len_str.parse().map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "Invalid CIGAR length")
                })?;
                mixed.push(MixedTracepointItem::CigarOp(len, op));
            } else {
                let parts: Vec<&str> = item.split(',').collect();
                if parts.len() == 2 {
                    let a: u64 = parts[0].parse().map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidData, "Invalid tracepoint value")
                    })?;
                    let b: u64 = parts[1].parse().map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidData, "Invalid tracepoint value")
                    })?;
                    mixed.push(MixedTracepointItem::Tracepoint(a, b));
                }
            }
        }
        Ok((TracepointData::Mixed(mixed), TracepointType::Mixed))
    } else {
        let has_single_values = items.iter().any(|s| !s.contains(','));
        if has_single_values {
            let mut variable = Vec::new();
            for item in items {
                if item.contains(',') {
                    let parts: Vec<&str> = item.split(',').collect();
                    if parts.len() == 2 {
                        let a: u64 = parts[0].parse().map_err(|_| {
                            io::Error::new(io::ErrorKind::InvalidData, "Invalid tracepoint value")
                        })?;
                        let b: u64 = parts[1].parse().map_err(|_| {
                            io::Error::new(io::ErrorKind::InvalidData, "Invalid tracepoint value")
                        })?;
                        variable.push((a, Some(b)));
                    }
                } else {
                    let a: u64 = item.parse().map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidData, "Invalid tracepoint value")
                    })?;
                    variable.push((a, None));
                }
            }
            Ok((TracepointData::Variable(variable), TracepointType::Variable))
        } else {
            let mut pairs = Vec::new();
            for item in items {
                let parts: Vec<&str> = item.split(',').collect();
                if parts.len() == 2 {
                    let a: u64 = parts[0].parse().map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidData, "Invalid tracepoint value")
                    })?;
                    let b: u64 = parts[1].parse().map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidData, "Invalid tracepoint value")
                    })?;
                    pairs.push((a, b));
                }
            }
            Ok((TracepointData::Standard(pairs), TracepointType::Standard))
        }
    }
}

// ============================================================================
// SEEKABLE READER WITH INDEX
// ============================================================================

/// Index for O(1) random access to records in a BPAF file
#[derive(Debug)]
pub struct BpafIndex {
    /// File offset for each record (byte position in .bpaf file)
    offsets: Vec<u64>,
}

impl BpafIndex {
    const INDEX_MAGIC: &'static [u8; 4] = b"BPAI";

    /// Save index to .bpaf.idx file
    pub fn save(&self, idx_path: &str) -> io::Result<()> {
        let mut file = File::create(idx_path)?;

        // Write magic and version
        file.write_all(Self::INDEX_MAGIC)?;
        file.write_all(&[1u8])?; // Version 1

        // Write number of offsets
        write_varint(&mut file, self.offsets.len() as u64)?;

        // Write all offsets
        for &offset in &self.offsets {
            write_varint(&mut file, offset)?;
        }

        Ok(())
    }

    /// Load index from .bpaf.idx file
    pub fn load(idx_path: &str) -> io::Result<Self> {
        let file = File::open(idx_path)?;
        let mut reader = BufReader::new(file);

        // Read and verify magic
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != Self::INDEX_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid index magic",
            ));
        }

        // Read version
        let mut version = [0u8; 1];
        reader.read_exact(&mut version)?;
        if version[0] != 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Unsupported index version: {}", version[0]),
            ));
        }

        // Read number of offsets
        let num_offsets = read_varint(&mut reader)? as usize;

        // Read all offsets
        let mut offsets = Vec::with_capacity(num_offsets);
        for _ in 0..num_offsets {
            offsets.push(read_varint(&mut reader)?);
        }

        Ok(Self { offsets })
    }

    /// Get number of records in index
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    /// Check if index is empty
    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }
}

/// Build index by scanning through BPAF file
pub fn build_index(bpaf_path: &str) -> io::Result<BpafIndex> {
    info!("Building index for {}", bpaf_path);

    let mut file = File::open(bpaf_path)?;
    let header = BinaryPafHeader::read(&mut file)?;

    let mut offsets = Vec::with_capacity(header.num_records as usize);

    // Skip codecs if adaptive format
    let strategy = CompressionStrategy::from_code(header.flags)?;
    if strategy.requires_codecs() {
        // Skip position codec
        AdaptiveCodec::read(&mut file)?;
        // Skip score codec
        AdaptiveCodec::read(&mut file)?;
    }

    // Record offset of each record
    for _ in 0..header.num_records {
        offsets.push(file.stream_position()?);
        skip_record(&mut file, strategy.requires_codecs())?;
    }

    info!("Index built: {} records", offsets.len());
    Ok(BpafIndex { offsets })
}

/// Skip a record without parsing (for building index)
fn skip_record<R: Read + Seek>(reader: &mut R, _is_adaptive: bool) -> io::Result<()> {
    // Skip core PAF fields (varints)
    read_varint(reader)?; // query_name_id
    read_varint(reader)?; // query_start
    read_varint(reader)?; // query_end
    reader.seek(SeekFrom::Current(1))?; // strand
    read_varint(reader)?; // target_name_id
    read_varint(reader)?; // target_start
    read_varint(reader)?; // target_end
    read_varint(reader)?; // residue_matches
    read_varint(reader)?; // alignment_block_len
    reader.seek(SeekFrom::Current(1))?; // mapping_quality
    reader.seek(SeekFrom::Current(1))?; // tp_type
    reader.seek(SeekFrom::Current(1))?; // complexity_metric
    read_varint(reader)?; // max_complexity

    // Skip tracepoints
    let num_items = read_varint(reader)? as usize;
    if num_items > 0 {
        // Skip compressed position data
        let pos_len = read_varint(reader)? as usize;
        reader.seek(SeekFrom::Current(pos_len as i64))?;

        // Skip compressed score data
        let score_len = read_varint(reader)? as usize;
        reader.seek(SeekFrom::Current(score_len as i64))?;
    }

    // Skip tags
    let num_tags = read_varint(reader)? as usize;
    for _ in 0..num_tags {
        reader.seek(SeekFrom::Current(2))?; // key
        let mut tag_type = [0u8; 1];
        reader.read_exact(&mut tag_type)?;
        match tag_type[0] {
            b'i' => reader.seek(SeekFrom::Current(8))?,
            b'f' => reader.seek(SeekFrom::Current(4))?,
            b'Z' => {
                let len = read_varint(reader)? as usize;
                reader.seek(SeekFrom::Current(len as i64))?
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid tag type",
                ))
            }
        };
    }

    Ok(())
}

/// Seekable reader for BPAF files with O(1) random access
pub struct BpafReader {
    file: File,
    index: BpafIndex,
    header: BinaryPafHeader,
    string_table: StringTable,
    // Adaptive format codecs (if applicable)
    first_val_codec: Option<AdaptiveCodec>,
    second_val_codec: Option<AdaptiveCodec>,
}

impl BpafReader {
    /// Open a BPAF file with index (builds index if .bpaf.idx doesn't exist)
    pub fn open(bpaf_path: &str) -> io::Result<Self> {
        let idx_path = format!("{}.idx", bpaf_path);

        // Load or build index
        let index = if Path::new(&idx_path).exists() {
            info!("Loading existing index: {}", idx_path);
            BpafIndex::load(&idx_path)?
        } else {
            info!("No index found, building...");
            let idx = build_index(bpaf_path)?;
            idx.save(&idx_path)?;
            debug!("Index saved to {}", idx_path);
            idx
        };

        // Open file and read header
        let mut file = File::open(bpaf_path)?;
        let header = BinaryPafHeader::read(&mut file)?;

        // Load codecs if strategy requires them
        let strategy = CompressionStrategy::from_code(header.flags)?;
        let (first_val_codec, second_val_codec) = if strategy.requires_codecs() {
            let pos = AdaptiveCodec::read(&mut file)?;
            let score = AdaptiveCodec::read(&mut file)?;
            (Some(pos), Some(score))
        } else {
            (None, None)
        };

        // Lazy-load string table only when needed (expensive for large files)
        let string_table = StringTable::new();

        Ok(Self {
            file,
            index,
            header,
            string_table,
            first_val_codec,
            second_val_codec,
        })
    }

    /// Open a BPAF file without index (for offset-based access only)
    ///
    /// Use this if you have your own offset storage (like impg) and only need:
    /// - get_alignment_record_at_offset()
    /// - get_tracepoints_at_offset()
    ///
    /// This skips index loading entirely - much faster open time.
    pub fn open_without_index(bpaf_path: &str) -> io::Result<Self> {
        // Open file and read header
        let mut file = File::open(bpaf_path)?;
        let header = BinaryPafHeader::read(&mut file)?;

        // Load codecs if strategy requires them
        let strategy = CompressionStrategy::from_code(header.flags)?;
        let (first_val_codec, second_val_codec) = if strategy.requires_codecs() {
            let pos = AdaptiveCodec::read(&mut file)?;
            let score = AdaptiveCodec::read(&mut file)?;
            (Some(pos), Some(score))
        } else {
            (None, None)
        };

        // Empty index - not used for offset-based access
        let index = BpafIndex {
            offsets: Vec::new(),
        };
        let string_table = StringTable::new();

        Ok(Self {
            file,
            index,
            header,
            string_table,
            first_val_codec,
            second_val_codec,
        })
    }

    /// Load string table (call this if you need sequence names)
    pub fn load_string_table(&mut self) -> io::Result<()> {
        if !self.string_table.is_empty() {
            return Ok(()); // Already loaded
        }

        // Seek to string table (at end of file, after all records)
        if !self.index.offsets.is_empty() {
            self.file.seek(SeekFrom::Start(
                self.index.offsets[self.index.offsets.len() - 1],
            ))?;
            skip_record(&mut self.file, self.header.flags & FLAG_ADAPTIVE != 0)?;
        }

        self.string_table = StringTable::read(&mut self.file)?;
        Ok(())
    }

    /// Get number of records
    pub fn len(&self) -> usize {
        self.index.len()
    }

    /// Check if reader is empty
    pub fn is_empty(&self) -> bool {
        self.index.is_empty()
    }

    /// Get header information
    pub fn header(&self) -> &BinaryPafHeader {
        &self.header
    }

    /// Get string table (loads on first access if needed)
    pub fn string_table(&mut self) -> io::Result<&StringTable> {
        if self.string_table.is_empty() {
            self.load_string_table()?;
        }
        Ok(&self.string_table)
    }

    /// Get immutable reference to string table (must be loaded first with load_string_table)
    pub fn string_table_ref(&self) -> &StringTable {
        &self.string_table
    }

    /// Get full alignment record by ID - O(1) random access
    pub fn get_alignment_record(&mut self, record_id: u64) -> io::Result<AlignmentRecord> {
        if record_id >= self.index.len() as u64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Record ID {} out of range (max: {})",
                    record_id,
                    self.index.len() - 1
                ),
            ));
        }

        let offset = self.index.offsets[record_id as usize];
        self.get_alignment_record_at_offset(offset)
    }

    /// Get alignment record by file offset (for impg compatibility)
    pub fn get_alignment_record_at_offset(&mut self, offset: u64) -> io::Result<AlignmentRecord> {
        self.file.seek(SeekFrom::Start(offset))?;

        let strategy = CompressionStrategy::from_code(self.header.flags)?;
        match strategy {
            CompressionStrategy::Varint => AlignmentRecord::read(&mut self.file),
            CompressionStrategy::VarintRaw => read_record_varint(&mut self.file, false),
            CompressionStrategy::Huffman => read_record_huffman(
                &mut self.file,
                self.first_val_codec.as_ref().unwrap(),
                self.second_val_codec.as_ref().unwrap(),
                true,
            ),
            CompressionStrategy::HuffmanRaw => read_record_huffman(
                &mut self.file,
                self.first_val_codec.as_ref().unwrap(),
                self.second_val_codec.as_ref().unwrap(),
                false,
            ),
            CompressionStrategy::VarintHuffman => read_record_varint_huffman(
                &mut self.file,
                self.first_val_codec.as_ref().unwrap(),
                self.second_val_codec.as_ref().unwrap(),
                true,
            ),
            CompressionStrategy::VarintHuffmanRaw => read_record_varint_huffman(
                &mut self.file,
                self.first_val_codec.as_ref().unwrap(),
                self.second_val_codec.as_ref().unwrap(),
                false,
            ),
            CompressionStrategy::Smart => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Smart strategy should not be stored in file",
            )),
        }
    }

    /// Get tracepoints only (optimized) - O(1) random access by record ID
    /// Returns: (tracepoints, tp_type, complexity_metric, max_complexity)
    ///
    /// Optimized for tracepoint-only access:
    /// - Default compression: skips unnecessary fields
    /// - Adaptive compression: reads full record but avoids field extraction overhead
    pub fn get_tracepoints(
        &mut self,
        record_id: u64,
    ) -> io::Result<(TracepointData, TracepointType, ComplexityMetric, u64)> {
        if record_id >= self.index.len() as u64 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "Record ID {} out of range (max: {})",
                    record_id,
                    self.index.len() - 1
                ),
            ));
        }

        let offset = self.index.offsets[record_id as usize];
        self.get_tracepoints_at_offset(offset)
    }

    /// Get tracepoints by file offset (for impg compatibility)
    /// Returns: (tracepoints, tp_type, complexity_metric, max_complexity)
    ///
    /// Use this if you have stored the actual file offsets - skips index lookup.
    pub fn get_tracepoints_at_offset(
        &mut self,
        offset: u64,
    ) -> io::Result<(TracepointData, TracepointType, ComplexityMetric, u64)> {
        self.file.seek(SeekFrom::Start(offset))?;

        let strategy = CompressionStrategy::from_code(self.header.flags)?;
        if strategy.requires_codecs() {
            // Codec-based strategies - must read full record
            let record = match strategy {
                CompressionStrategy::Huffman => read_record_huffman(
                    &mut self.file,
                    self.first_val_codec.as_ref().unwrap(),
                    self.second_val_codec.as_ref().unwrap(),
                    true,
                ),
                CompressionStrategy::HuffmanRaw => read_record_huffman(
                    &mut self.file,
                    self.first_val_codec.as_ref().unwrap(),
                    self.second_val_codec.as_ref().unwrap(),
                    false,
                ),
                CompressionStrategy::VarintHuffman => read_record_varint_huffman(
                    &mut self.file,
                    self.first_val_codec.as_ref().unwrap(),
                    self.second_val_codec.as_ref().unwrap(),
                    true,
                ),
                CompressionStrategy::VarintHuffmanRaw => read_record_varint_huffman(
                    &mut self.file,
                    self.first_val_codec.as_ref().unwrap(),
                    self.second_val_codec.as_ref().unwrap(),
                    false,
                ),
                _ => unreachable!("Strategy should require codecs"),
            }?;
            Ok((
                record.tracepoints,
                record.tp_type,
                record.complexity_metric,
                record.max_complexity,
            ))
        } else {
            // Varint-only strategies - can skip fields
            read_varint(&mut self.file)?; // query_name_id - SKIP
            read_varint(&mut self.file)?; // query_start - SKIP
            read_varint(&mut self.file)?; // query_end - SKIP
            self.file.seek(SeekFrom::Current(1))?; // strand - SKIP
            read_varint(&mut self.file)?; // target_name_id - SKIP
            read_varint(&mut self.file)?; // target_start - SKIP
            read_varint(&mut self.file)?; // target_end - SKIP
            read_varint(&mut self.file)?; // residue_matches - SKIP
            read_varint(&mut self.file)?; // alignment_block_len - SKIP
            self.file.seek(SeekFrom::Current(1))?; // mapping_quality - SKIP

            // Read only what we need
            let mut tp_type_buf = [0u8; 1];
            self.file.read_exact(&mut tp_type_buf)?;
            let tp_type = TracepointType::from_u8(tp_type_buf[0])
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            let mut metric_buf = [0u8; 1];
            self.file.read_exact(&mut metric_buf)?;
            let complexity_metric = complexity_metric_from_u8(metric_buf[0])?;

            let max_complexity = read_varint(&mut self.file)?;

            // Read tracepoints
            let tracepoints = AlignmentRecord::read_tracepoints(&mut self.file, tp_type)?;

            Ok((tracepoints, tp_type, complexity_metric, max_complexity))
        }
    }

    /// Iterator over all records (sequential access)
    pub fn iter_records(&mut self) -> RecordIterator<'_> {
        RecordIterator {
            reader: self,
            current_id: 0,
        }
    }
}

/// Iterator for sequential record access
pub struct RecordIterator<'a> {
    reader: &'a mut BpafReader,
    current_id: u64,
}

impl<'a> Iterator for RecordIterator<'a> {
    type Item = io::Result<AlignmentRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current_id >= self.reader.len() as u64 {
            return None;
        }

        let result = self.reader.get_alignment_record(self.current_id);
        self.current_id += 1;
        Some(result)
    }
}
