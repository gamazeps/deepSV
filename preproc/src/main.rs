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
use std::fs;
use std::io::BufRead;
use std::io::BufReader;
use std::sync::{Arc, Mutex};

use consts::FULL_1000GP_VCF_PATH;
use vcf_record::{parse_vcf_file};
use generate_read::{generate_reads};

use std::time::{Instant};

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
        let l = line.expect("should be able to read a line");
        let _ = fs::create_dir(format!("../data/supporting_reads/{}", l.clone()));
        whitelist.insert(l);
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

    let mut size = 0;
    for record in records {
        sender.send(record).unwrap();
        size += 1;
    }
    let size = size;

    let n_threads = 32;
    let mut thread_ids = Vec::with_capacity(n_threads);

    let beg = Instant::now();
    for i in 0..n_threads {
        let recv = receiver.clone();
        let id = thread::spawn(move || {
            let mut cnt = 0;
            loop {
                let record = recv.lock().unwrap().try_recv();
                match record {
                    Ok(r) => {
                        //println!("Start: thread {}, record {}, time {:?}",
                        //         i, sample, Instant::now());
                        generate_reads(r);
                        //println!("End: thread {}, record {}, time {:?}",
                        //         i, sample, Instant::now());
                        cnt+=1;
                        if (cnt % 20) == 0 {
                            println!("thread {} processed {} variants of the total {}, it is {:.3}%",
                                     i, cnt, size, ((n_threads * cnt) as f64) / (size as f64));
                        }
                    }
                    Err(err) => {
                        println!("{:?}", err);
                        break;
                    }
                }
            }
        });
        thread_ids.push(id);
    }

    for id in thread_ids {
        id.join().expect("should not have failed");
    }
    let end = Instant::now();
    println!("{:?}", end.duration_since(beg));
}
