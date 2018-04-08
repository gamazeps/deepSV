"""
Let's write this shit as if it were flume !
"""
import collections
import json
import glob
import logging
import os
import sys

import pysam

import utils

median_insert_size = 400
median_read_size = 101
window = median_insert_size + median_read_size


def record_dict(record):
    # Needed because record.info is a VariantRecordInfo object which
    # overrides its __getitem__ method in order to act as a dict.
    # Since the VariantRecordInfo class is not serializable, we need to do it ourselve.
    return {
       "chrom": record.chrom,
       "start": record.start,
       "stop": record.stop,
       "id": record.id,
       "alts": record.alts,
       "info": dict(record.info.items())
    }


def extract_reads(records, sample, conf, bam_path=None):
    """
    The `bam_path` argument is here as a hack to allow me to work locally
    """
    if bam_path is None:
        bam_dir = os.path.join(conf["alignments"], sample, "alignment")
        if not os.path.exists(bam_dir):
            logging.info("No reads for {}".format(sample))
            return
        fnames = glob.glob(os.path.join(bam_dir, "*.bam"))
        assert(len(fnames) == 1)
        bamfile = pysam.AlignmentFile(fnames[0], 'rb')
    else:
        bamfile = pysam.AlignmentFile(bam_path, 'rb')

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
            f.write(json.dumps(record_dict(record)))
    logging.info("Finished extracting reads for {}".format(sample))


def main():
    utils.set_logging()
    logging.info("Starting to extract reads")

    if len(sys.argv) != 2:
        print("Please provide 1 argument: configuration.json")
        sys.exit(1)

    conf = utils.get_json(sys.argv[1])

    logging.info("Starting to parse the VCF")
    vcf = pysam.VariantFile(conf['vcf_path'], 'rb')

    per_sample = collections.defaultdict(list)
    for record in vcf.fetch('1'):
        # We want to remove variants with uncertain locations
        if "CIEND" in record.info.keys() or "CIPOS" in record.info.keys():
            continue
        per_sample[record.info["SAMPLE"]].append(record)
    logging.info("Finished parsing the VCF")

    #na12878 = per_sample["NA12878"]
    for (sample, records) in per_sample.items():
        extract_reads(records, sample, conf)

    logging.info("Finished extracting the reads")


if __name__ == "__main__":
    main()
