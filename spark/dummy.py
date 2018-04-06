"""
Let's write this shit as if it were flume !
"""
import collections
import json
import matplotlib.pyplot as plt

import pysam

median_insert_size = 400
median_read_size = 101
window = median_insert_size + median_read_size

def plot_histogram(records):
    per_sample = collections.defaultdict(int)
    for record in records:
        per_sample[record.info["SAMPLE"]] += 1
    x = [v for (k, v) in per_sample.items()]
    plt.hist(x, bins=100)
    plt.show()


def record_dict(record):
    return {
       "chrom": record.chrom,
       "start": record.start,
       "stop": record.stop,
       "id": record.id,
       "alts": record.alts,
       "info": dict(record.info.items())
    }


def extract_reads(records):
    fname = "../data/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
    samfile = pysam.AlignmentFile(fname, 'rb')
    for record in records[:1]:
        variant_size = record.stop - record.stop
        supporting_reads = pysam.AlignmentFile("{}.bam".format(record.id),
                                               "wb",
                                               template=samfile)

        # Taking two halves, but probably not needed at that step, could be used only at
        # the tensor building phase.
        # Could even build the tensor there...
        # NOTE: the sam file was poorly generated once, but I cannot reproduce.
        if variant_size > (2 * window + median_read_size + 20):
            l_reads = samfile.fetch(record.chrom, record.start - window, record.start + window)
            r_reads = samfile.fetch(record.chrom, record.stop - window, record.stop + window)
            for read in l_reads:
                supporting_reads.write(read)
            for read in r_reads:
                supporting_reads.write(read)
        else:
            for read in samfile.fetch(record.chrom, record.start - window, record.stop + window):
                supporting_reads.write(read)
        supporting_reads.close()

        fafile = pysam.FastaFile("human_g1k_v37.fasta")
        with open("{}.fa".format(record.id), 'w') as f:
            f.write("# {} {} {}\n".format(record.chrom,
                                          record.start - window,
                                          record.start + window))
            f.write("{}\n".format(str(fafile.fetch(record.chrom,
                                                   record.start - window,
                                                   record.start + window))))
            f.write("# {} {} {}\n".format(record.chrom,
                                          record.stop - window,
                                          record.stop + window))
            f.write("{}\n".format(str(fafile.fetch(record.chrom,
                                                   record.stop - window,
                                                   record.stop + window))))

        with open("{}.json".format(record.id), 'w') as f:
            f.write(json.dumps(record_dict(record)))


def main():
    fname = "estd219.GRCh37.variant_call.vcf.gz"
    vcf = pysam.VariantFile(fname, 'rb')

    per_sample = collections.defaultdict(list)
    for record in vcf.fetch('1'):
        # We want to remove variants with uncertain locations
        if "CIEND" in record.info.keys() or "CIPOS" in record.info.keys():
            continue
        per_sample[record.info["SAMPLE"]].append(record)

    na12878 = per_sample["NA12878"]
    print(len(na12878))
    extract_reads(na12878)

if __name__ == "__main__":
    main()
