import glob
import json
import logging
import multiprocessing
import pysam
import random
import sys
import time
import os
import numpy as np
import tensorflow as tf
import argparse

from read_pairs import ReadPair, RefSeq
import utils
import deepsv_tensor

def find_variant_files(path, sample):
    fnames = glob.glob("{}/{}/*bam".format(path, sample))
    return [fname[:-4] for fname in fnames]


def build_tensor(basename, reads_limit=150, ref_rep=5):
    metadata = utils.get_json(basename + ".json")
    read_pairs = build_read_pairs(basename + ".bam", reads_limit - ref_rep)
    begin = time.time()
    ref = RefSeq(basename + ".fa")
    encoder = deepsv_tensor.TensorEncoder(n_channels=4)

    tensor = deepsv_tensor.DeepSVTensor(encoder=encoder,
                                        metadata=metadata,
                                        pairs_capacity=reads_limit)

    for _ in range(5):
        tensor.insert_read_pair(ref)

    for pair in read_pairs:
        tensor.insert_read_pair(pair)

    return tensor


def build_read_pairs(fname, sample_limit):
    with pysam.AlignmentFile(fname, 'rb') as bamfile:
        records = dict()
        for read in bamfile.fetch(until_eof=True):
            if not read.is_paired:
                continue
            if read.qname not in records:
                records[read.qname] = ReadPair(read.qname)
            records[read.qname].add_read(read)

        # Here we sample up to `sample_limit` reads at random
        records_list = [v for (_, v) in records.items()]
        n_records = len(records_list)
        sampled_indices = sorted(random.sample(range(n_records), min(n_records, sample_limit)))
        sampled_records_list = [records_list[i] for i in sampled_indices]

        return sorted(sampled_records_list, key=lambda rp: rp.leftmost())


def get_whitelist(fname):
    with open(fname, "r") as f:
        return [l.strip() for l in f]


def process_sample(sample, conf):
    names = find_variant_files(conf["supporting_reads_path"], sample)
    if len(names) == 0:
        return

    logging.info("Starting processing {}".format(sample))

    times = np.zeros((len(names),), dtype='f')

    out_path = os.path.join(conf["tensors_path"], sample + '.tfrecords')

    # Open a TFRecordWriter for the output-file.
    with tf.python_io.TFRecordWriter(out_path) as writer:

        for i, fname in enumerate(names):
            start = time.time()
            tensor = build_tensor(fname)

            # Create a dict with the data we want to save in the
            # TFRecords file. You can add more relevant data here.
            # Wrap the data as TensorFlow Features.
            data = \
                {
                    'data': tf.train.Feature(bytes_list=tf.train.BytesList(value=[tensor.tensor.tostring()])),
                    'label': tf.train.Feature(int64_list=tf.train.Int64List(value=[tensor.label()]))
                }

            feature = tf.train.Features(feature=data)

            # Wrap again as a TensorFlow Example.
            example = tf.train.Example(features=feature)

            # Serialize the data.
            serialized = example.SerializeToString()

            # Write the serialized data to the TFRecords file.
            writer.write(serialized)
            times[i] = time.time() - start

    logging.info("done saving {} to {}".format(sample, out_path))
    logging.info('avg time: {} per tensor'.format(times.mean()))

    # The return value is needed for multiprocessing
    return None


def main():
    parser = argparse.ArgumentParser(description='Tool for generating TFRecords from bam files')
    parser.add_argument("--conf", type=str, required=True)
    parser.add_argument("--whitelist", type=str, required=True)
    parser.add_argument("--log-output", type=str, default=None)
    args = parser.parse_args()

    utils.set_logging(fname=args.log_output)

    conf_fname = args.conf
    whitelist_fname = args.whitelist

    conf = utils.get_json(conf_fname)
    samples = get_whitelist(whitelist_fname)

    if not os.path.exists(conf['tensors_path']):
        logging.info("Creating the {} directory".format(conf['tensors_path']))
        os.makedirs(conf['tensors_path'])

    samples_size = len(samples)
    logging.info("There is a total of {} samples to process".format(samples_size))

    n_threads = conf.get("n_threads", 1)

    if n_threads > 1:
        p = multiprocessing.Pool(n_threads)
        p.starmap(process_sample, ((sample, conf) for sample in samples))
    else:
        for sample in samples:
            process_sample(sample, conf)

    logging.info("Finished generating tensors")


if __name__ == "__main__":
    main()
