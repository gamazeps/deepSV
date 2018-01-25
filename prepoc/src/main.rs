/// File standards follow the description of https://github.com/samtools/hts-specs
/// according to commit: 36576b330c81b0f40c37c97c6aa4529b488e48a1

mod picard_stats;
mod vcf_record;

use std::io;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;

use vcf_record::{VCFRecord, parse_record};

fn read_file() -> Result<(), io::Error> {
    let f = try!(File::open("../data/PHASE3_SV_NA12878.vcf"));
    //let f = try!(File::open("/data/fraimund/ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV/vcf/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf"));
    let file = BufReader::new(&f);

    let records: Vec<VCFRecord> = file.lines()
                                      .filter_map(|line| parse_record(line.unwrap()))
                                      .collect();

    let count = records.len();
    println!("Records: {}", count);

    let ci_count = records.iter().filter(|record| record.has_ci()).count();
    println!("Inconfident Records: {}", ci_count);
    println!("Confident Records: {}", count - ci_count);

    Ok(())
}

fn main() {
    let _ = read_file();
}
