/// File standards follow the description of https://github.com/samtools/hts-specs
/// according to commit: 36576b330c81b0f40c37c97c6aa4529b488e48a1

use std::collections::HashMap;
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
    info: HashMap<String, InfoField>
}

/// VCF record according to the V
impl VCFRecord {
    fn new(chromosome: String, pos: Option<u64>, id: String, reference: String,
           alt: Vec<VariantType>, quality: Option<u32>, filter: String,
           info: HashMap<String, InfoField>)
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

    fn is_alt(&self, alts: Vec<VariantType>) -> bool {
        self.alt == alts
    }
}

/// Fields here follow the ones from the
/// [released dbvar file](ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV/vcf/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf.gz)
#[derive(Debug, PartialEq, Eq)]
enum VariantType {
    /// Copy number variable region.
    CNV,
    /// Deletion relative to the reference.
    DEL,
    /// Region of elevated copy number relative to the reference.
    DUP, 
    /// Insertion of sequence relative to the reference.
    INS, 
    /// Inversion of reference sequence.
    INV, 
    /// Undocumented.
    SVA, 
    /// Undocumented.
    ALU,
    /// Undocumented.
    ME,
    /// Undocumented.
    L1,
    /// Human Endogenous RetroVirus.
    HERV
}

fn parse_alt(alt: &str) -> Vec<VariantType> {
    let pattern: &[_] = &['<', '>'];
    let fields = alt.trim_matches(pattern).split(":");
    fields.filter_map(|v| match v {
        "CNV"  => Some(VariantType::CNV),
        "DEL"  => Some(VariantType::DEL),
        "DUP"  => Some(VariantType::DUP),
        "INS"  => Some(VariantType::INS),
        "INV"  => Some(VariantType::INV),
        "SVA"  => Some(VariantType::SVA),
        "ALU"  => Some(VariantType::ALU),
        "ME"   => Some(VariantType::ME),
        "L1"   => Some(VariantType::L1),
        "HERV" => Some(VariantType::HERV),
        other  => {println!("Unknown ALT value: {}", other); None}
    }).collect()
}

/// Fields here follow the ones from the
/// [released dbvar file](ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV/vcf/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf.gz)
#[derive(Debug, PartialEq, Eq)]
enum InfoField {
    /// ID is a dbVar accession.
    DBVARID, 

    /// Confidence interval around END for imprecise variants.
    CIEND(i64, u64), 

    /// Confidence interval around POS for imprecise variants.
    CIPOS(i64, u64),

    /// Any additional information about this call (free text, enclose in double quotes).
    DESC(String),

    /// End position of the variant described in this record.
    END(u64),  

    /// Imprecise structural variation.
    IMPRECISE, 

    /// For imprecise variants, if END represents an inner_stop or outer_stop coordinate, use a
    /// comma-delimited list to indicate this; e.g., if END is an inner_stop and its value is
    /// 108442336, you would enter '108442336,.'.
    ENDrange(u64, u64),

    /// For imprecise variants, if POS represents an inner_start or outer_start coordinate, use a
    /// comma-delimited list to indicate this; e.g., if POS is an inner_start and its value is
    /// 2865734, you would enter '.,2865734'.
    POSrange(u64, u64),

    /// Type of structural variant, must be one of: DEL, INS, DUP, INV, CNV.
    SVTYPE(String),

    /// Submitted variant call id.
    CALLID(String),

    /// The parent variant region accession(s).
    REGION(String),

    /// The experiment_id (from EXPERIMENTS tab) of the experiment that was used to generate this
    /// call.
    EXPERIMENT(u64),

    /// Sample_id from dbVar submission; every call must have SAMPLE or SAMPLESET, but NOT BOTH.
    SAMPLE(String)
}

fn parse_info(info: &str) -> HashMap<String, InfoField> {
    let fields = info.split(";");
    fields.fold(HashMap::new(), |mut m, v| {
        let pattern: &[_] = &['=', ','];
        let infos : Vec<&str> = v.split(pattern).collect();
        match infos[0] {
        "DBVARID"    => m.insert(infos[0].to_string(),
                                 InfoField::DBVARID),
        "CIEND"      => m.insert(infos[0].to_string(),
                                 InfoField::CIEND(i64::from_str(infos[1]).unwrap(),
                                                  u64::from_str(infos[2]).unwrap())),
        "CIPOS"      => m.insert(infos[0].to_string(),
                                 InfoField::CIPOS(i64::from_str(infos[1]).unwrap(),
                                                  u64::from_str(infos[2]).unwrap())),
        "DESC"       => m.insert(infos[0].to_string(),
                                 InfoField::DESC(String::from_str(infos[1]).unwrap())),
        "END"        => m.insert(infos[0].to_string(),
                                 InfoField::END(u64::from_str(infos[1]).unwrap())),
        "IMPRECISE"  => m.insert(infos[0].to_string(), InfoField::IMPRECISE),
        "ENDrange"   => m.insert(infos[0].to_string(), 
                                 InfoField::ENDrange(u64::from_str(infos[1]).unwrap(),
                                                     u64::from_str(infos[2]).unwrap())),
        "POSrange"   => m.insert(infos[0].to_string(),
                                 InfoField::POSrange(u64::from_str(infos[1]).unwrap(),
                                                     u64::from_str(infos[2]).unwrap())),
        "SVTYPE"     => m.insert(infos[0].to_string(),
                                 InfoField::SVTYPE(String::from_str(infos[1]).unwrap())),
        "CALLID"     => m.insert(infos[0].to_string(),
                                 InfoField::CALLID(String::from_str(infos[1]).unwrap())),
        "REGION"     => m.insert(infos[0].to_string(),
                                 InfoField::REGION(String::from_str(infos[1]).unwrap())),
        "EXPERIMENT" => m.insert(infos[0].to_string(),
                                 InfoField::EXPERIMENT(u64::from_str(infos[1]).unwrap())),
        "SAMPLE"     => m.insert(infos[0].to_string(),
                                        InfoField::SAMPLE(String::from_str(infos[1]).unwrap())),
        other        => {println!("Unknown INFO field: {}", other); None}
        };
        m
    })
}

fn parse_record(input: String) -> Option<VCFRecord> {
    // We remove headers that we do not parse.
    if input.starts_with("#") {
        return None
    }

    let words : Vec<&str> = input.split('\t').collect();
    Some(VCFRecord::new(
        String::from_str(words[0]).unwrap(),
        Some(u64::from_str(words[1]).unwrap()),
        String::from_str(words[2]).unwrap(),
        String::from_str(words[3]).unwrap(),
        parse_alt(words[4]),
        u32::from_str(words[5]).ok(),
        String::from_str(words[6]).unwrap(),
        parse_info(words[7]),
    ))
}

fn read_file() -> Result<(), io::Error> {
    let f = try!(File::open("../data/PHASE3_SV_NA12878.vcf"));
    //let f = try!(File::open("/data/fraimund/ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV/vcf/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf"));
    let file = BufReader::new(&f);

    let records = file.lines().filter_map(|line| parse_record(line.unwrap()));

    let records = records.filter(|record| record.is_alt(vec![VariantType::DEL]));
    
    for r in records.take(20) {
        println!("{:?}", r);
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
