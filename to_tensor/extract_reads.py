import collections
import json
import glob
import logging
import multiprocessing
import os
import sys
import argparse

import pysam

import utils

median_insert_size = 400
median_read_size = 101
window = median_insert_size + median_read_size


class VariantRecordWrapper(object):
    # Needed because record.info is a VariantRecordInfo object which
    # overrides its __getitem__ method in order to act as a dict.
    # Since the VariantRecordInfo class is not serializable, we need to do it ourselve.
    def __init__(self, record):
        self.chrom = record.chrom
        self.start = record.start
        self.stop = record.stop
        self.id = record.id
        self.alts = record.alts
        self.info = dict(record.info.items())

    def to_json(self):
        return {
            'chrom': self.chrom,
            'start': self.start,
            'stop': self.stop,
            'id': self.id,
            'alts': self.alts,
            'info': self.info
        }


def extract_reads(sample, records, conf):
    bam_dir = os.path.join(conf["alignments"], sample, "alignment")
    if not os.path.exists(bam_dir):
        return
    fnames = glob.glob(os.path.join(bam_dir, "*.bam"))
    assert(len(fnames) == 1)
    bamfile = pysam.AlignmentFile(fnames[0], 'rb')

    ref_fa = pysam.FastaFile(conf["reference_path"])
    out_dir = os.path.join(conf['supporting_reads_path'], sample)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for record in records:
        bam_fname = os.path.join(out_dir, "{}.bam".format(record.id))
        supporting_reads = pysam.AlignmentFile(bam_fname, "wb", template=bamfile)

        variant_size = record.stop - record.stop
        # Taking two halves, but probably not needed at that step, could be used only at
        # the tensor building phase.
        # Could even build the tensor there...
        # NOTE: the sam file was poorly generated once, but I cannot reproduce.
        if variant_size > (2 * window + median_read_size + 20):
            l_reads = bamfile.fetch(record.chrom, record.start - window, record.start + window)
            r_reads = bamfile.fetch(record.chrom, record.stop - window, record.stop + window)
            for read in l_reads:
                supporting_reads.write(read)
            for read in r_reads:
                supporting_reads.write(read)
        else:
            for read in bamfile.fetch(record.chrom, record.start - window, record.stop + window):
                supporting_reads.write(read)
        supporting_reads.close()


        fa_fname = os.path.join(out_dir, "{}.fa".format(record.id))
        with open(fa_fname, 'w') as f:
            f.write("# {} {} {}\n".format(record.chrom,
                                          record.start - window,
                                          record.start + window))
            f.write("{}\n".format(str(ref_fa.fetch(record.chrom,
                                                   record.start - window,
                                                   record.start + window))))
            f.write("# {} {} {}\n".format(record.chrom,
                                          record.stop - window,
                                          record.stop + window))
            f.write("{}\n".format(str(ref_fa.fetch(record.chrom,
                                                   record.stop - window,
                                                   record.stop + window))))

        json_fname = os.path.join(out_dir, "{}.json".format(record.id))
        with open(json_fname, 'w') as f:
            f.write(json.dumps(record.to_json()))
    logging.info("Finished extracting reads for {}".format(sample))
    return None


def main():
    parser = argparse.ArgumentParser(description='Tool for extracting reads from bam files')
    parser.add_argument("--conf", type=str, required=True)
    args = parser.parse_args()
    conf = utils.get_json(args.conf)

    utils.set_logging()
    logging.info("Starting to extract reads")

    logging.info("Starting to parse the VCF")
    vcf = pysam.VariantFile(conf['vcf_path'], 'rb')

    per_sample = collections.defaultdict(list)
    for record in vcf.fetch():
        # We want to remove variants with uncertain locations
        if "CIEND" in record.info.keys() or "CIPOS" in record.info.keys():
            continue
        per_sample[record.info["SAMPLE"]].append(VariantRecordWrapper(record))
    logging.info("Finished parsing the VCF")

    n_threads = conf.get("n_threads", 1)
    logging.info("Starting to generate the reads on {} threads".format(n_threads))

    if n_threads > 1:
        p = multiprocessing.Pool(n_threads)
        p.starmap(extract_reads,
                  ((sample, records, conf) for (sample, records) in per_sample.items()))
    else:
        for (sample, records) in per_sample.items():
            extract_reads(sample, records, conf)

    logging.info("Finished extracting the reads")


if __name__ == "__main__":
    main()
