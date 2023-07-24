use vcf::{VCFReader, VCFWriter, VCFError, U8Vec};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufReader, BufRead, stdout};
use clap::Parser;
use sha1::{Sha1, Digest};

#[derive(Parser)]
#[command(name = "VCF ID adder")]
#[command(author = "Erik Garrison <erik.garrison@gmail.com>")]
#[command(version = "0.1")]
#[command(about = "Sets IDs on records using the position and alleles to define the ID.", long_about = None)]
struct Cli {
    #[arg(short,long)]
    input: String,
    #[arg(short,long)]
    sha1_hash: bool, // use a hash value instead of the position and alleles
    #[arg(short,long,default_value = "_")]
    delim: String, // delimiter for the ID field
}

fn main() -> Result<(), VCFError> {

    let cli = Cli::parse();

    let use_hash_id = cli.sha1_hash;

    let delim = cli.delim;

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

    while reader.next_record(&mut vcf_record)? {
        let chrom = std::str::from_utf8(&vcf_record.chromosome)?;
        let pos = vcf_record.position.to_string();
        // concatenate the alternative alleles into a string, separated by ,
        let alternative = vcf_record.alternative.iter().map(|x| std::str::from_utf8(x)).collect::<Result<Vec<&str>, _>>()?.join(",");
        let reference = std::str::from_utf8(&vcf_record.reference)?.to_string();
        let mut new_id = format!("{}{}{}{}{}{}{}", chrom, delim, pos, delim, reference, delim, alternative);
        if use_hash_id {
            let mut hasher = Sha1::new();
            hasher.update(new_id.as_bytes());
            // take id as a string that's the first 4 bytes of the hash
            new_id = hasher.finalize().iter().take(4).map(|b| format!("{:02x}", b)).collect::<String>();
        }
        vcf_record.id = vec![U8Vec::from(new_id.as_bytes())];
        writer.write_record(&vcf_record)?;
    }

    Ok(())
}

