use std::io;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;

struct VCFRecord {
    chromosome: String,
    pos: u64,
    id: String,
    reference: String,
    alt: String,
    quality: u32,
    filter: String,
    info: String
}


fn read_file() -> Result<(), io::Error> {
    let f = try!(File::open("../data/PHASE3_SV_NA12878.vcf"));
    let file = BufReader::new(&f);

    let mut i : u32 = 0;

    for line in file.lines() {
        i += 1;
        let l = line.unwrap();
        let words : Vec<&str> = l.split('\t').collect();
        println!("{:?}", words); 
        if i > 10 {
            break;
        }
    }
    Ok(())
}

fn main() {
    let _ = read_file();
}
