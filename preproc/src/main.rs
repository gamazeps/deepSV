/// File standards follow the description of https://github.com/samtools/hts-specs
/// according to commit: 36576b330c81b0f40c37c97c6aa4529b488e48a1

#[macro_use]
extern crate serde_derive;
extern crate glob;

mod picard_stats;
mod vcf_record;
mod generate_read;
mod config;

use std::env;
use std::collections::HashSet;
use std::sync::mpsc::channel;
use std::thread;
use std::fs::File;
use std::fs;
use std::io::BufRead;
use std::io::BufReader;
use std::sync::{Arc, Mutex};

use vcf_record::{parse_vcf_file};
use generate_read::{generate_reads};
use config::Config;

use std::time::{Instant};

const USE_STR : &'static str= "Use of preproc: cargo run CONFIG TARGETS";

fn get_config() -> Config {
    let args: Vec<String> = env::args().collect();
    let fname = match args.len() {
        3 => args[1].clone(),
        _ => panic!(USE_STR),
    };
    Config::from_name(fname)
}

fn whitelisted_samples(destination: String) -> HashSet<String> {
    let args: Vec<String> = env::args().collect();
    let fname = match args.len() {
        3 => args[2].clone(),
        _ => panic!(USE_STR),
    };

    let f = File::open(fname).expect("Unable to open the given whitelist file");
    let file = BufReader::new(&f);

    let mut whitelist = HashSet::new();
    for line in file.lines().into_iter() {
        let l = line.expect("should be able to read a line");
        let _ = fs::create_dir(format!("{}/{}", destination, l.clone()));
        whitelist.insert(l);
    }

    whitelist
}

fn main() {
    let config = get_config();

    let records = parse_vcf_file(config.vcf_path.clone());
    let whitelist = whitelisted_samples(config.destination_dir.clone());

    println!("Opening the VCF file");
    let records = records
        .into_iter()
        .filter(move |record| record.is_simple()
                && !record.has_ci()
                && whitelist.contains(&record.sample()));

    let (sender, receiver) = channel();
    let receiver = Arc::new(Mutex::new(receiver));

    println!("Sending reads to the channels");
    let mut size = 0;
    for record in records {
        sender.send(record).unwrap();
        size += 1;
    }
    let size = size;

    let n_threads = 32;
    let mut thread_ids = Vec::with_capacity(n_threads);

    println!("Launching the threads");
    let beg = Instant::now();
    for i in 0..config.n_threads {
        let recv = receiver.clone();
        let conf = config.clone();
        let id = thread::spawn(move || {
            let mut cnt = 0;
            loop {
                let record = recv.lock().unwrap().try_recv();
                match record {
                    Ok(r) => {
                        generate_reads(r, &conf);
                        cnt+=1;
                        if (cnt % 20) == 0 {
                            println!("thread {} processed {} variants of the total {},\
                                     it is {:.3}% in {}s",
                                     i, cnt, size,
                                     100.0 * ((n_threads * cnt) as f64) / (size as f64),
                                     Instant::now().duration_since(beg).as_secs());
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
