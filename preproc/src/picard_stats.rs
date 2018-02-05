use std::process::Command;
use std::thread;
use std::thread::JoinHandle;
use std::sync::Arc;

#[allow(non_snake_case, dead_code)]
pub struct CollectAlignmentSummaryMetrics {
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

#[allow(non_snake_case, dead_code)]
pub struct InsertSizeMetrics {
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

#[allow(dead_code)]
pub struct PicardMetrics {
    insert_size_metrics: InsertSizeMetrics,
    alignement_metrics: CollectAlignmentSummaryMetrics
}

impl PicardMetrics {
    fn new(insert_size_metrics: InsertSizeMetrics,
           alignement_metrics: CollectAlignmentSummaryMetrics)
        -> PicardMetrics {
        PicardMetrics {
            alignement_metrics: alignement_metrics,
            insert_size_metrics: insert_size_metrics
        }
    }
}

fn insert_metrics(picard_path: Arc<String>, bam_file: Arc<String>, reference_file: Arc<String>)
    -> JoinHandle<Option<InsertSizeMetrics>> {

    let insert_metrics_file = format!("{}.insert.tmp", bam_file);
    thread::spawn(move || {
        let _ = Command::new("java")
            .arg("-jar")
            .arg(picard_path.as_ref())
            .arg("InsertSizeMetrics")
            .arg(format!("R={}", reference_file))
            .arg(format!("I={}", bam_file))
            .arg(format!("O={}", insert_metrics_file))
            .arg("H=/tmp/who_cares.pdf")
            .arg("M=0.5")
            .spawn();
        None
    })
}

fn alignement_metrics(picard_path: Arc<String>, bam_file: Arc<String>, reference_file: Arc<String>)
    -> JoinHandle<Option<CollectAlignmentSummaryMetrics>> {

    let alignement_metrics_file = format!("{}.allignment.tmp", bam_file);
    thread::spawn(move || {
        let _ = Command::new("java")
            .arg("-jar")
            .arg(picard_path.as_ref())
            .arg("CollectAlignmentSummaryMetrics")
            .arg(format!("R={}", reference_file))
            .arg(format!("I={}", bam_file))
            .arg(format!("O={}", alignement_metrics_file));
        None
    })
}

#[allow(dead_code)]
pub fn picard_metrics(picard_path: String, bam_file: String, reference_file: String)
    -> PicardMetrics {

    let picard_path = Arc::new(picard_path);
    let bam_file = Arc::new(bam_file);
    let reference_file = Arc::new(reference_file);

    let alignement_handle = alignement_metrics(picard_path.clone(),
                                               bam_file.clone(),
                                               reference_file.clone());
    let insert_handle = insert_metrics(picard_path.clone(),
                                       bam_file.clone(),
                                       reference_file.clone());

    let insert = insert_handle.join().unwrap() // Result from join.
                                     .unwrap(); // Option of the function.
    let alignment = alignement_handle.join().unwrap() // Option from join.
                                            .unwrap(); // Option from the function.

    PicardMetrics::new(insert, alignment)
}

