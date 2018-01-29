/// File standards follow the description of https://github.com/samtools/hts-specs
/// according to commit: 36576b330c81b0f40c37c97c6aa4529b488e48a1

mod picard_stats;
mod vcf_record;

use vcf_record::{group_by_sample, parse_vcf_file};

fn main() {
    //let records = parse_vcf_file("../data/PHASE3_SV_NA12878.vcf".to_owned());
    let records = parse_vcf_file(
        "/data/fraimund/ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study\
        /estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV/vcf\
        /estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37\
        .submitted.variant_call.germline.vcf".to_owned()
    );

    let count = records.len();
    println!("Records: {}", count);

    for (k, v) in group_by_sample(records).into_iter() {
        println!("{}: {} variants", k, v.len());
    }

    //let ci_count = records.iter().filter(|record| record.has_ci()).count();
    //println!("Inconfident Records: {}", ci_count);
    //println!("Confident Records: {}", count - ci_count);
}
