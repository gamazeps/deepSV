/// Module for calling the nice function that will generat the reads using samtools
/// in a nice file.

use std::io::prelude::*;
use std::fs::File;
use std::process::Command;

use glob::glob;

use vcf_record::{InfoField, VCFRecord};
use consts::{REFERENCE_FA};

extern crate serde;
extern crate serde_json;

// TODO(gamazeps) do not hardcode the naming of the files for the samples
pub fn generate_reads(record: VCFRecord) {
    let mut fnames : Vec<_> = glob("../data/alignment/*.bam")
        .expect("Failed to read glob pattern")
        .collect();
    if fnames.len() == 0 {
        return ();
    }
    if fnames.len() > 1 {
        println!("There are two files, {:?}", fnames);
        return ();
    }

    let fname = fnames.pop().unwrap().unwrap();

    let median_insert_size = 400;
    let median_read_size = 100;
    let window = median_insert_size + median_read_size;

    if !record.has_ci() {
        let end = match record.get_info("END".to_owned()) {
            Some(InfoField::END(e)) => e,
            _ => panic!("END field should contain a pos")
        };
        let pos = record.pos().unwrap();

        // TODO(gamazeps): this is a horrible interraction with the borrowchecker.
        let mut c: Command = Command::new("samtools");
        c.arg("view")
         .arg(format!("{}", fname.display()))
         .arg("-M")
         .arg(format!("{}:{}-{}",
                      record.chromosome(), pos - window, pos + window))
         .arg(format!("{}:{}-{}",
                      record.chromosome(), end - window, end + window));

        let output = c.output().expect("Failed to execute samtools view");

        let sample = record.sample();
        let sv_type = match record.get_info("SVTYPE".to_owned()) {
            Some(InfoField::SVTYPE(e)) => e,
            _ => panic!("SVTYPE field should contain a type")
        };

        let fname_sam = fname_ext(sample.clone(), record.id(), sv_type.clone(),
                                  record.pos().unwrap(), end, "sam");

        // TODO(gamazeps): use https://github.com/rust-lang/rust/pull/42133/files for the output
        // This is currently shady as fuck...
        let mut buffer = File::create(fname_sam.clone()).expect("Failed to create {}");
        buffer.write(&output.stdout).expect("Failed to write the data to {}");
        buffer.sync_all().expect("Failed to sync {} to disk");

        let mut ref_c: Command = Command::new("samtools");
        ref_c.arg("faidx")
         .arg(REFERENCE_FA)
         .arg(format!("{}:{}-{}",
                      record.chromosome(), pos - window, pos + window))
         .arg(format!("{}:{}-{}",
                      record.chromosome(), end - window, end + window));
        let output = ref_c.output().expect("Failed to execute samtools faidx");
        let fname_fa = fname_ext(sample.clone(), record.id(), sv_type.clone(),
                                 record.pos().unwrap(), end, "fa");
        let mut buffer = File::create(fname_fa.clone()).expect("Failed to create {}");
        buffer.write(&output.stdout).expect("Failed to write the data to {}");
        buffer.sync_all().expect("Failed to sync {} to disk");

        let fname_json = fname_ext(sample, record.id(), sv_type,
                                   record.pos().unwrap(), end, "json");
        let metadata = VarianrtMetadata::new(record, fname_fa, fname_sam);
        metadata.write_to_file(fname_json);
    }
}

fn fname_ext(sample: String, id: String, sv_type: String, pos: u64, end: u64, ext: &str) -> String {
        format!(
            "../data/supporting_reads/{}/{}.{}.{}-{}.{}",
            sample,
            id,
            sv_type,
            pos,
            end,
            ext
        )
}

#[derive(Serialize, Deserialize)]
struct VarianrtMetadata {
    record: VCFRecord,
    fa_location: String,
    sam_location: String
}

impl VarianrtMetadata {
    fn new(record: VCFRecord, fa_location: String, sam_location: String) -> VarianrtMetadata {
        VarianrtMetadata {
            record: record,
            fa_location: fa_location,
            sam_location: sam_location
        }
    }

    fn write_to_file(&self, fname: String) {
        let json = serde_json::to_string(self).expect("unable to serialize metadata");
        let mut buffer = File::create(fname).expect("Failed to create metadata file");
        buffer.write(json.as_bytes()).expect("Failed to write the metadata");
        buffer.sync_all().expect("Failed to sync metadata to disk");
    }
}
