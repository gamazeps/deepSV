/// Module for calling the nice function that will generat the reads using samtools
/// in a nice file.

use vcf_record::{InfoField, VCFRecord};

use std::process::Command;

// TODO(gamazeps) do not hardcode the naming of the files for the samples
pub fn generate_reads_for_na12878(record: VCFRecord) -> Option<Command> {

    let na12878_bam_path = "NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam";
    // TODO(gamazeps): use https://github.com/rust-lang/rust/pull/42133/files for the output

    assert_eq!(record.get_info("SAMPLE".to_owned()), Some(InfoField::SAMPLE("NA12878".to_owned())));

    if record.has_ci() {
        None
    } else {
        let end = match record.get_info("END".to_owned()) {
            Some(InfoField::END(e)) => e,
            _ => panic!("END fieled should contain a pos")
        };

        // TODO(gamazeps): this is a horrible interraction with the borrowchecker.
        let mut c: Command = Command::new("samtools");
        c.arg("view")
         .arg(na12878_bam_path)
         .arg(format!("{}:{}-{}", record.chromosome(), record.pos().unwrap(), end));
        println!("{:?}", c);
        Some(c)
    }


    //samtools view  "1:4204668-4204717" > toto

}
