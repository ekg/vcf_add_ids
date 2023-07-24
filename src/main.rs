use vcf::{VCFReader, VCFWriter, VCFError, U8Vec};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufReader, BufRead, stdout, Write};
use clap::Parser;
use sha1::{Sha1, Digest};

#[derive(Parser)]
#[command(name = "VCF ID adder")]
#[command(author = "Erik Garrison <erik.garrison@gmail.com>")]
#[command(version = "0.1")]
#[command(about = "Sets IDs on records using the position and alleles to define the ID.", long_about = None)]
struct Cli {
    /// input VCF file
    #[arg(short,long)]
    input: String,
    /// take the first 4 bytes of a sha1 hash of the record pos/ref/alt as the id
    #[arg(short,long)]
    sha1_hash: bool,
    /// use this delimiter
    #[arg(short,long,default_value = "_")]
    delim: String, // delimiter for the ID field
    /// a prefix for the ids
    #[arg(short,long,default_value = "")]
    prefix: String,
}

struct Config {
    input: String,
    sha1_hash: bool,
    delim: String,
    prefix: String,
    output: Box<dyn Write>,
}

fn main() -> Result<(), VCFError> {
    let cli = Cli::parse();
    let config = Config {
        input: cli.input,
        sha1_hash: cli.sha1_hash,
        delim: cli.delim,
        prefix: cli.prefix,
        output: Box::new(stdout()),
    };
    modify_ids(config)
}

fn modify_ids(config: Config) -> Result<(), VCFError> {
    let filename = config.input;
    let file = File::open(filename.clone())?;
    let buf_reader: Box<dyn BufRead> = if filename.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut reader = VCFReader::new(buf_reader)?;
    let header = reader.header().clone();
    let mut writer = VCFWriter::new(config.output, &header)?;

    // prepare VCFRecord object
    let mut vcf_record = reader.empty_record();

    while reader.next_record(&mut vcf_record)? {
        let chrom = std::str::from_utf8(&vcf_record.chromosome)?;
        let pos = vcf_record.position.to_string();
        // concatenate the alternative alleles into a string, separated by ,
        let alternative = vcf_record.alternative.iter().map(|x| std::str::from_utf8(x)).collect::<Result<Vec<&str>, _>>()?.join(",");
        let reference = std::str::from_utf8(&vcf_record.reference)?.to_string();
        let mut new_id = format!("{}{}{}{}{}{}{}", chrom, config.delim, pos, config.delim, reference, config.delim, alternative);
        if config.sha1_hash {
            let mut hasher = Sha1::new();
            hasher.update(new_id.as_bytes());
            // take id as a string that's the first 4 bytes of the hash
            new_id = hasher.finalize().iter().take(4).map(|b| format!("{:02x}", b)).collect::<String>();
        }
        new_id = format!("{}{}", config.prefix, new_id);
        vcf_record.id = vec![U8Vec::from(new_id.as_bytes())];
        writer.write_record(&vcf_record)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_modify_ids() {
        let config = Config {
            input: "./test/z.vcf.gz".to_string(),
            sha1_hash: false,
            delim: "_".to_string(),
            prefix: "".to_string(),
            // write to /dev/null
            output: Box::new(File::create("/dev/null").unwrap()),
        };
        let result = modify_ids(config);
        assert!(result.is_ok());
    }
}
