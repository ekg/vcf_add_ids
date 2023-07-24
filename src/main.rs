use vcf::{VCFReader, VCFWriter, VCFError, VCFHeader, VCFRecord, U8Vec};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufReader, BufRead, stdout, Write};
use clap::Parser;

#[derive(Parser)]
#[command(name = "VCF ID adder")]
#[command(author = "Erik Garrison <erik.garrison@gmail.com>")]
#[command(version = "0.1")]
#[command(about = "Sets IDs on records using the position and alleles to define the ID.", long_about = None)]
struct Cli {
    #[arg(short,long)]
    input: String,
}

fn main() -> Result<(), VCFError> {

    let cli = Cli::parse();

    // copy cli.input into filenam
    let filename = cli.input;
    let file = File::open(filename.clone())?;
    let buf_reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut reader = VCFReader::new(buf_reader)?;
    let header = reader.header().clone();
    let mut writer = VCFWriter::new(stdout(), &header)?;

    // prepare VCFRecord object
    let mut vcf_record = reader.empty_record();

    while reader.next_record(&mut vcf_record).is_ok() {
        let chrom = std::str::from_utf8(&vcf_record.chromosome)?;
        let pos = vcf_record.position.to_string();
        // as a string
        let alternative = std::str::from_utf8(vcf_record.alternative.first().unwrap())?.to_string();
        let reference = std::str::from_utf8(&vcf_record.reference)?.to_string();
        let new_id = format!("{}:{}:{}:{}", chrom, pos, reference, alternative);
        vcf_record.id = vec![U8Vec::from(new_id.as_bytes())];
        writer.write_record(&vcf_record)?;
    }

    Ok(())
}

