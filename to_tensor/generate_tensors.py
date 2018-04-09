import cPickle
import glob
import h5py
import json
import logging
import multiprocessing
import pysam
import random
import sys
import os

from read_pairs import ReadPair, RefSeq
import utils
import deepsv_tensor

global_conf = None
count = 0

def find_variant_files(path, sample):
    fnames = glob.glob("{}/{}/*bam".format(path, sample))
    return [fname[:-4] for fname in fnames]


def build_tensor(basename, reads_limit=150):
    metadata = utils.get_json(basename + ".json")
    read_pairs = build_read_pairs(basename + ".bam", reads_limit)
    ref = RefSeq(basename + ".fa")
    encoder = deepsv_tensor.TensorEncoder(n_channels=8, sam_channels=4, ref_channels=4)

    tensor = deepsv_tensor.DeepSVTensor(encoder=encoder,
                                        metadata=metadata,
                                        pairs_capacity=reads_limit)
    tensor.insert_ref(ref)
    for pair in read_pairs:
        tensor.insert_read_pair(pair)

    return tensor


def build_read_pairs(fname, sample_limit):
    bamfile = pysam.AlignmentFile(fname, 'rb')
    records = dict()
    for read in bamfile.fetch(until_eof=True):
        if not read.is_paired:
            continue
        if read.qname not in records:
            records[read.qname] = ReadPair(read.qname)
        records[read.qname].add_read(read)

    # Here we sample up to `sample_limit` reads at random
    records_list = [v for (_, v) in records.iteritems()]
    n_records = len(records_list)
    sampled_indices = sorted(random.sample(range(n_records), min(n_records, sample_limit)))
    sampled_records_list = [records_list[i] for i in sampled_indices]

    return sorted(sampled_records_list, key=lambda rp: rp.leftmost())


def get_whitelist(fname):
    with open(fname, "r") as f:
        return [l.strip() for l in f]


def process_variant(fname, draw=False):
    tensor = build_tensor(fname)
    if draw:
        tensor.dummy_image("{}.png".format(fname), draw_bp=True)
        logging.info("{}.png".format(fname))
    return tensor


def process_sample(conf, sample):
    names = find_variant_files(conf["reads_path"], sample)
    # TODO(gamazeps): keeping all the tensors in memory uses a lot of RAM for nothing.
    # The tensors could be written one by one, however we have enough memory on the cluster
    # to do that with 32 threads in parrallel (750GB used at the max).
    # It could easily be optimized by writing the tensors one by one to hdf5.
    tensors = [process_variant(fname) for fname in names]
    logging.info("done processing {}".format(sample))

    if len(tensors) == 0:
        return

    # Needed for encoding the json metadata
    dt = h5py.special_dtype(vlen=bytes)

    raw_data = [t.tensor for t in tensors]
    labels = [t.label() for t in tensors]
    metadata = [json.dumps(t.metadata) for t in tensors]

    with h5py.File("{}/{}.hdf5".format(conf["tensors_path"], sample), "w") as f:
        data_dset= f.create_dataset("data", data=raw_data, compression="lzf", dtype="u1")
        labels_dset= f.create_dataset("labels", data=labels, compression="lzf", dtype="u1")
        metadata_dset= f.create_dataset("metadata", data=metadata, compression="lzf", dtype=dt)
    logging.info("done saving {} to hdf5".format(sample))


def par_process_sample(sample):
    process_sample(global_conf, sample)
    return 0


def main():
    # The global variables are ugly, but we need them to be global so that they
    # can de passed to multi_processing
    global global_conf

    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Please provide at least 2 arguments: configuration.json whitelist [log]")
        sys.exit(1)

    if len(sys.argv) == 4:
        print("logging to {}".format(sys.argv[3]))
        utils.set_logging(fname=sys.argv[3])
    else:
        utils.set_logging()

    conf_fname = sys.argv[1]
    whitelist_fname = sys.argv[2]

    global_conf = utils.get_json(conf_fname)
    samples = get_whitelist(whitelist_fname)

    if not os.path.exists(global_conf['tensors_path']):
        logging.info("Creating the {} directory".format(global_conf['tensors_path']))
        os.makedirs(global_conf['tensors_path'])

    samples_size = len(samples)
    logging.info("There is a total of {} samples to process".format(samples_size))

    n_threads = global_conf.get("n_threads", 1)

    if n_threads > 1:
        p = multiprocessing.Pool(n_threads)
        p.map(par_process_sample, samples)
    else:
        for sample in samples:
            process_sample(global_conf, sample)

    logging.info("Finished generating tensors")
    sys.exit(0)


if __name__ == "__main__":
    main()
