extern crate serde;
extern crate serde_json;

use std::fs::File;
use std::io::BufReader;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Config {
    pub destination_dir: String,
    pub data_dir: String,
    pub vcf_path: String,
    pub refence_path: String,
    pub n_threads: u64,
}

impl Config {
    pub fn from_name(fname: String) -> Config {
        let file = File::open(fname.clone())
            .expect(&format!("Failed to open the configuration file {}", fname.clone()));
        let buffer = BufReader::new(file);
        serde_json::from_reader(buffer)
            .expect(&format!("Failed to read or parse the configuration file {}", fname))
    }
}
