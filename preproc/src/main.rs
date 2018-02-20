/// File standards follow the description of https://github.com/samtools/hts-specs
/// according to commit: 36576b330c81b0f40c37c97c6aa4529b488e48a1

#[macro_use]
extern crate serde_derive;

mod picard_stats;
mod vcf_record;
mod generate_read;
mod consts;

use std::env;
use std::collections::HashSet;
use std::sync::mpsc::channel;
use std::thread;


use consts::NA12878_VCF_PATH;
use vcf_record::{parse_vcf_file};
use generate_read::{generate_reads_for_na12878};


fn whitelisted_samples() -> HashSet<String> {
    let args: Vec<String> = env::args().collect();
    let fname = match args.len() {
        2 => Some(args[1]),
        _ => panic!("Please provide one argument with a path to the whitelisted samples"),
    }

    let f = File::open(fname).expect("Unable to open the given whitelist file");
    let file = BufReader::new(&f);

    let mut whitelist = HashSet::new();
    for line in file.lines().into_iter() {
        whitelist.insert(line);
    }

    whitelist
}

fn main() {
    let records = parse_vcf_file(NA12878_VCF_PATH.to_owned());
    let whitelist = whitelisted_samples();

    let records = records
        .into_iter()
        .filter(move |record| record.is_simple()
                && !record.has_ci()
                && whitelist.contains(record));

    let (sender, receiver) = channel();

    for record in records {
        sender.send(record).unwrap();
    }

    let mut thread_ids = Vector::new();
    for i in 0..30 {
    }

    for id in thread_ids {
        id.join(),expect("should not have failed");
    }
}
