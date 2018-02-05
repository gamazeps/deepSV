/// Module for calling the nice function that will generat the reads using samtools
/// in a nice file.

use std::io::prelude::*;
use std::fs::File;
use std::process::Command;

use vcf_record::{InfoField, VCFRecord};
use consts::NA12878_BAM_PATH;


// TODO(gamazeps) do not hardcode the naming of the files for the samples
pub fn generate_reads_for_na12878(record: VCFRecord) {
    assert_eq!(
        record.get_info("SAMPLE".to_owned()),
        Some(InfoField::SAMPLE("NA12878".to_owned()))
    );

    if !record.has_ci() {
        let end = match record.get_info("END".to_owned()) {
            Some(InfoField::END(e)) => e,
            _ => panic!("END fieled should contain a pos")
        };

        // TODO(gamazeps): this is a horrible interraction with the borrowchecker.
        let mut c: Command = Command::new("samtools");
        c.arg("view")
         .arg(NA12878_BAM_PATH)
         .arg(format!("{}:{}-{}", record.chromosome(), record.pos().unwrap(), end));
        println!("{:?}", c);


        let output = c.output().expect("Failed to execute process");
        println!("status: {}", output.status);
        println!("stderr: {}", String::from_utf8_lossy(&output.stderr));

        // TODO(gamazeps): use https://github.com/rust-lang/rust/pull/42133/files for the output
        // This is currently shady as fuck...
        let mut buffer = File::create(format!("{}.sam", record.id()))
            .expect("should be able to create a file");
        buffer.write(&output.stdout);
        buffer.sync_all().expect("should be able to sync the data");
    }
}
