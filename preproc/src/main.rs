/// File standards follow the description of https://github.com/samtools/hts-specs
/// according to commit: 36576b330c81b0f40c37c97c6aa4529b488e48a1

#[macro_use]
extern crate serde_derive;
extern crate glob;

mod picard_stats;
mod vcf_record;
mod generate_read;
mod consts;

use std::env;
use std::collections::HashSet;
use std::sync::mpsc::channel;
use std::thread;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::sync::{Arc, Mutex};

use consts::FULL_1000GP_VCF_PATH;
use vcf_record::{parse_vcf_file};
use generate_read::{generate_reads};


fn whitelisted_samples() -> HashSet<String> {
    let args: Vec<String> = env::args().collect();
    let fname = match args.len() {
        2 => args[1].clone(),
        _ => panic!("Please provide one argument with a path to the whitelisted samples"),
    };

    let f = File::open(fname).expect("Unable to open the given whitelist file");
    let file = BufReader::new(&f);

    let mut whitelist = HashSet::new();
    for line in file.lines().into_iter() {
        whitelist.insert(line.expect("should be able to read a line"));
    }

    whitelist
}

fn main() {
    let records = parse_vcf_file(FULL_1000GP_VCF_PATH.to_owned());
    let whitelist = whitelisted_samples();

    let records = records
        .into_iter()
        .filter(move |record| record.is_simple()
                && !record.has_ci()
                && whitelist.contains(&record.sample()));

    let (sender, receiver) = channel();
    let receiver = Arc::new(Mutex::new(receiver));

    for record in records {
        sender.send(record).unwrap();
    }

    let mut thread_ids = Vec::new();

    for i in 0..30 {
        let recv = receiver.clone();
        let id = thread::spawn(move || {
            let mut cnt = 0;
            loop {
                let record = recv.lock().unwrap().recv();
                match record {
                    Ok(r) => generate_reads(r),
                    Err(err) => {
                        println!("{:?}", err);
                        break;
                    }
                }
                cnt+=1;
                if (cnt % 100) == 0 {
                    println!("Thread {} processed {} reads", i, cnt);
                }
            }
        });
        thread_ids.push(id);
    }

    for id in thread_ids {
        id.join().expect("should not have failed");
    }
}
