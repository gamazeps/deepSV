/// File standards follow the description of https://github.com/samtools/hts-specs
/// according to commit: 36576b330c81b0f40c37c97c6aa4529b488e48a1

use std::collections::HashMap;
use std::io;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::str::FromStr;

// TODO split that shit

struct CollectAlignmentSummaryMetrics {
    /// One of either UNPAIRED (for a fragment run), FIRST_OF_PAIR when metrics are for only the
    /// first read in a paired run, SECOND_OF_PAIR when the metrics are for only the second read
    /// in a paired run or PAIR when the metrics are aggregated for both first and second reads
    /// in a pair.
    CATEGORY : u32,

    /// The total number of reads including all PF and non-PF reads. When CATEGORY equals PAIR this
    /// value will be 2x the number of clusters.
    TOTAL_READS : u32,

    /// The number of PF reads where PF is defined as passing Illumina's filter.
    PF_READS : u32,

    /// The fraction of reads that are PF (PF_READS / TOTAL_READS)
    PCT_PF_READS : u32,

    /// The number of PF reads that are marked as noise reads. A noise read is one which is
    /// composed entirely of A bases and/or N bases. These reads are marked as they are usually
    /// artifactual and are of no use in downstream analysis.
    PF_NOISE_READS : u32,

    /// The number of PF reads that were aligned to the reference sequence. This includes reads
    /// that aligned with low quality (i.e. their alignments are ambiguous).
    PF_READS_ALIGNED : u32,

    /// The percentage of PF reads that aligned to the reference sequence. PF_READS_ALIGNED /
    /// PF_READS
    PCT_PF_READS_ALIGNED : u32,

    /// The total number of aligned bases, in all mapped PF reads, that are aligned to the
    /// reference sequence.
    PF_ALIGNED_BASES : u32,

    /// The number of PF reads that were aligned to the reference sequence with a mapping quality of
    /// Q20 or higher signifying that the aligner estimates a 1/100 (or smaller) chance that the
    /// alignment is wrong.
    PF_HQ_ALIGNED_READS : u32,

    /// The number of bases aligned to the reference sequence in reads that were mapped at high
    /// quality. Will usually approximate PF_HQ_ALIGNED_READS * READ_LENGTH but may differ when
    /// either mixed read lengths are present or many reads are aligned with gaps.
    PF_HQ_ALIGNED_BASES : u32,

    /// The subset of PF_HQ_ALIGNED_BASES where the base call quality was Q20 or higher.
    PF_HQ_ALIGNED_Q20_BASES : u32,

    /// The median number of mismatches versus the reference sequence in reads that were aligned to
    /// the reference at high quality (i.e. PF_HQ_ALIGNED READS).
    PF_HQ_MEDIAN_MISMATCHES : u32,

    /// The rate of bases mismatching the reference for all bases aligned to the reference
    /// sequence.
    PF_MISMATCH_RATE : u32,

    /// The fraction of bases that mismatch the reference in PF HQ aligned reads.
    PF_HQ_ERROR_RATE : u32,

    /// The number of insertion and deletion events per 100 aligned bases. Uses the number of
    /// events as the numerator, not the number of inserted or deleted bases.
    PF_INDEL_RATE : u32,

    /// The mean read length of the set of reads examined. When looking at the data for a single
    /// lane with equal length reads this number is just the read length. When looking at data
    /// for merged lanes with differing read lengths this is the mean read length of all reads.
    MEAN_READ_LENGTH : u32,

    /// The number of aligned reads whose mate pair was also aligned to the reference.
    READS_ALIGNED_IN_PAIRS : u32,

    /// The fraction of reads whose mate pair was also aligned to the reference.
    /// READS_ALIGNED_IN_PAIRS / PF_READS_ALIGNED
    PCT_READS_ALIGNED_IN_PAIRS : u32,

    /// The number of (primary) aligned reads that are **not** "properly" aligned in pairs (as per
    /// SAM flag 0x2).
    PF_READS_IMPROPER_PAIRS : u32,

    /// The fraction of (primary) reads that are *not* "properly" aligned in pairs (as per SAM flag
    /// 0x2). PF_READS_IMPROPER_PAIRS / PF_READS_ALIGNED
    PCT_PF_READS_IMPROPER_PAIRS : u32,

    /// The number of instrument cycles in which 80% or more of base calls were no-calls.
    BAD_CYCLES : u32,

    /// The number of PF reads aligned to the positive strand of the genome divided by the number
    /// of PF reads aligned to the genome.
    STRAND_BALANCE : u32,

    /// The fraction of reads that map outside of a maximum insert size (usually 100kb) or that
    /// have the two ends mapping to different chromosomes.
    PCT_CHIMERAS : u32,

    /// The fraction of PF reads that are unaligned and match to a known adapter sequence right from
    /// the start of the read.
    PCT_ADAPTER : u32,
}

struct InsertSizeMetrics {
    /// The MEDIAN insert size of all paired end reads where both ends mapped to the same
    /// chromosome.
    MEDIAN_INSERT_SIZE : u32,

    /// The median absolute deviation of the distribution. If the distribution is essentially
    /// normal then the standard deviation can be estimated as ~1.4826 * MAD.
    MEDIAN_ABSOLUTE_DEVIATION : u32,

    /// The minimum measured insert size. This is usually 1 and not very useful as it is likely
    /// artifactual.
    MIN_INSERT_SIZE : u32,

    /// The maximum measure insert size by alignment. This is usually very high representing either
    /// an artifact or possibly the presence of a structural re-arrangement.
    MAX_INSERT_SIZE : u32,

    /// The mean insert size of the "core" of the distribution. Artefactual outliers in the
    /// distribution often cause calculation of nonsensical mean and stdev values. To avoid this
    /// the distribution is first trimmed to a "core" distribution of +/- N median absolute
    /// deviations around the median insert size. By default N=10, but this is configurable.
    MEAN_INSERT_SIZE : u32,

    /// Standard deviation of insert sizes over the "core" of the distribution.
    STANDARD_DEVIATION : u32,

    /// The total number of read pairs that were examined in the entire distribution.
    READ_PAIRS : u32,

    /// The pair orientation of the reads in this data category.
    PAIR_ORIENTATION : u32,

    /// The "width" of the bins, centered around the median, that encompass 10% of all read pairs.
    WIDTH_OF_10_PERCENT : u32,

    /// The "width" of the bins, centered around the median, that encompass 20% of all read pairs.
    WIDTH_OF_20_PERCENT : u32,

    /// The "width" of the bins, centered around the median, that encompass 30% of all read pairs.
    WIDTH_OF_30_PERCENT : u32,

    /// The "width" of the bins, centered around the median, that encompass 40% of all read pairs.
    WIDTH_OF_40_PERCENT : u32,

    /// The "width" of the bins, centered around the median, that encompass 50% of all read pairs.
    WIDTH_OF_50_PERCENT : u32,

    /// The "width" of the bins, centered around the median, that encompass 60% of all read pairs.
    WIDTH_OF_60_PERCENT : u32,

    /// The "width" of the bins, centered around the median, that encompass 70% of all read pairs.
    /// This metric divided by 2 should approximate the standard deviation when the insert size
    /// distribution is a normal distribution.
    WIDTH_OF_70_PERCENT : u32,

    /// The "width" of the bins, centered around the median, that encompass 80% of all read pairs.
    WIDTH_OF_80_PERCENT : u32,

    /// The "width" of the bins, centered around the median, that encompass 90% of all read pairs.
    WIDTH_OF_90_PERCENT : u32,

    /// The "width" of the bins, centered around the median, that encompass 100% of all read pairs.
    WIDTH_OF_99_PERCENT : u32,
}


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

    fn has_ci(&self) -> bool {
        self.info.contains_key("CIEND") || self.info.contains_key("CIPOS")
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
    //let f = try!(File::open("../data/PHASE3_SV_NA12878.vcf"));
    let f = try!(File::open("/data/fraimund/ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV/vcf/estd219_1000_Genomes_Consortium_Phase_3_Integrated_SV.GRCh37.submitted.variant_call.germline.vcf"));
    let file = BufReader::new(&f);

    let records: Vec<VCFRecord> = file.lines().filter_map(|line| parse_record(line.unwrap())).collect();

    let count = records.len();
    println!("Records: {}", count);

    let ci_count = records.iter().filter(|record| record.has_ci()).count();
    println!("Inconfident Records: {}", ci_count);
    println!("Confident Records: {}", count - ci_count);

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
