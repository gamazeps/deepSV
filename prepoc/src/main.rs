/// File standards follow the description of https://github.com/samtools/hts-specs
/// according to commit: 36576b330c81b0f40c37c97c6aa4529b488e48a1

mod picard_stats;
mod vcf_record;
mod generate_read;
mod consts;

use consts::NA12878_VCF_PATH;
use vcf_record::{parse_vcf_file};
use generate_read::{generate_reads_for_na12878};

fn main() {
    let records = parse_vcf_file(NA12878_VCF_PATH.to_owned());

    let count = records.len();
    println!("Records: {}", count);

    for r in records.into_iter().take(3) {
        generate_reads_for_na12878(r);
    }

    //let ci_count = records.iter().filter(|record| record.has_ci()).count();
    //println!("Inconfident Records: {}", ci_count);
    //println!("Confident Records: {}", count - ci_count);
}
