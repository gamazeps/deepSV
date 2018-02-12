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

    for r in records.into_iter().filter(|record| record.is_simple() && !record.has_ci()).take(20) {
        generate_reads_for_na12878(r);
    }
}
