use std::io;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::str::FromStr;

#[derive(Debug, PartialEq, Eq)]
struct VCFRecord {
    chromosome: String,
    pos: Option<u64>,
    id: String,
    reference: String,
    alt: Vec<VariantType>,
    quality: Option<u32>,
    filter: String,
    info: Vec<InfoField>
}

impl VCFRecord {
    fn new(chromosome: String, pos: Option<u64>, id: String, reference: String,
           alt: Vec<VariantType>, quality: Option<u32>, filter: String, info: Vec<InfoField>)
        -> VCFRecord {
        VCFRecord {
            chromosome: chromosome,
            pos: pos,
            id: id,
            reference: reference,
            alt: alt,
            quality: quality,
            filter: filter,
            info: info
        }
    }
}

#[derive(Debug)]
enum VariantType {
    CNV,
    DEL,
    DUP,
    INS,
    INV,
    SVA,
    ALU,
    ME,
    L1,
}

fn parse_alt(alt: &str) -> Vec<VariantType> {
    let pattern: &[_] = &['<', '>'];
    let fields = alt.trim_matches(pattern).split(":");
    fields.filter_map(|v| match v {
        "CNV" => Some(VariantType::CNV),
        "DEL" => Some(VariantType::DEL),
        "DUP" => Some(VariantType::DUP),
        "INS" => Some(VariantType::INS),
        "INV" => Some(VariantType::INV),
        "SVA" => Some(VariantType::SVA),
        "ALU" => Some(VariantType::ALU),
        "ME"  => Some(VariantType::ME),
        "L1"  => Some(VariantType::L1),
        other => {println!("Unknown ALT value: {}", other); None}
    }).collect()
}

/// Info fields from the dbvar database
#[derive(Debug)]
enum InfoField {
    DBVARID,
    CIEND(i64, u64),
    CIPOS(i64, u64),
    DESC(String),
    END(u64),
    IMPRECISE,
    ENDrange(u64, u64),
    POSrange(u64, u64),
    SVTYPE(String),
    CALLID(String),
    REGION(String),
    EXPERIMENT(u64),
    SAMPLE(String)
}

fn parse_info(info: &str) -> Vec<InfoField> {
    let fields = info.split(";");
    fields.filter_map(|v| {
        let pattern: &[_] = &['=', ','];
        let infos : Vec<&str> = v.split(pattern).collect();
        match infos[0] {
        "DBVARID"    => Some(InfoField::DBVARID),
        "CIEND"      => Some(InfoField::CIEND(i64::from_str(infos[1]).unwrap(),
                                              u64::from_str(infos[2]).unwrap())),
        "CIPOS"      => Some(InfoField::CIPOS(i64::from_str(infos[1]).unwrap(),
                                              u64::from_str(infos[2]).unwrap())),
        "DESC"       => Some(InfoField::DESC(String::from_str(infos[1]).unwrap())),
        "END"        => Some(InfoField::END(u64::from_str(infos[1]).unwrap())),
        "IMPRECISE"  => Some(InfoField::IMPRECISE),
        "ENDrange"   => Some(InfoField::ENDrange(u64::from_str(infos[1]).unwrap(),
                                                 u64::from_str(infos[2]).unwrap())),
        "POSrange"   => Some(InfoField::POSrange(u64::from_str(infos[1]).unwrap(),
                                                 u64::from_str(infos[2]).unwrap())),
        "SVTYPE"     => Some(InfoField::SVTYPE(String::from_str(infos[1]).unwrap())),
        "CALLID"     => Some(InfoField::CALLID(String::from_str(infos[1]).unwrap())),
        "REGION"     => Some(InfoField::REGION(String::from_str(infos[1]).unwrap())),
        "EXPERIMENT" => Some(InfoField::EXPERIMENT(u64::from_str(infos[1]).unwrap())),
        "SAMPLE"     => Some(InfoField::SAMPLE(String::from_str(infos[1]).unwrap())),
        other        => {println!("Unknown INFO field: {}", other); None}
        }
    }).collect()
}


fn parse_record(input: String) -> VCFRecord {
    let words : Vec<&str> = input.split('\t').collect();
    VCFRecord::new(
        String::from_str(words[0]).unwrap(),
        Some(u64::from_str(words[1]).unwrap()),
        String::from_str(words[2]).unwrap(),
        String::from_str(words[3]).unwrap(),
        parse_alt(words[4]),
        u32::from_str(words[5]).ok(),
        String::from_str(words[6]).unwrap(),
        parse_info(words[7]),
    )
}

fn read_file() -> Result<(), io::Error> {
    let f = try!(File::open("../data/PHASE3_SV_NA12878.vcf"));
    let file = BufReader::new(&f);

    for line in file.lines() {
        let l = line.unwrap();
        let _ = parse_record(l);
    }
    Ok(())
}

fn main() {
    let _ = read_file();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_info() {
        assert_eq!(
            vec![InfoField::DBVARID,
                 InfoField::CALLID("DEL_pindel_91_NA12878".to_string()),
                 InfoField::SVTYPE("DEL".to_string()),
                 InfoField::EXPERIMENT(9),
                 InfoField::SAMPLE("NA12878".to_string()),
                 InfoField::END(2911850),
                 InfoField::REGION("esv3818169".to_string())],
            parse_info("DBVARID;CALLID=DEL_pindel_91_NA12878;SVTYPE=DEL;EXPERIMENT=9;SAMPLE=NA12878;END=2911850;REGION=esv3818169"))
    }
}
