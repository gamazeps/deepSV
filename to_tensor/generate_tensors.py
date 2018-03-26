import glob
import h5py
from multiprocessing import Pool
import cPickle
import sys
import random
import logging
import json

from read_pairs import SamRead, ReadPair, RefSeq
import utils
import deepsv_tensor
import disk_storage

global_conf = None
count = 0

def find_variant_files(path, sample):
    fnames = glob.glob("{}/{}/*sam".format(path, sample))
    return [fname[:-4] for fname in fnames]


def build_tensor(basename, reads_limit=150):
    metadata = utils.get_json(basename + ".json")
    record = metadata["record"]
    read_pairs = build_read_pairs(basename + ".sam", reads_limit)
    ref = RefSeq(basename + ".fa")
    encoder = deepsv_tensor.TensorEncoder(n_channels=8, sam_channels=4, ref_channels=4)

    tensor = deepsv_tensor.DeepSVTensor(encoder=encoder, metadata=record, pairs_capacity=reads_limit)
    tensor.insert_ref(ref)
    for pair in read_pairs:
        tensor.insert_read_pair(pair)

    return tensor


def build_read_pairs(fname, sample_limit):
    with open(fname, "r") as f:
        records = dict()
        for record in f:
            read = SamRead(record)
            if not (read.flag & int("0xc0", base=16)):
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
    tensors = [process_variant(fname) for fname in names]
    logging.info("done processing {}".format(sample))

    # Needed for encoding the json metadata
    dt = h5py.special_dtype(vlen=bytes)

    raw_data = [t.tensor for t in tensors]
    labels = [t.label() for t in tensors]
    metadata = [json.dumps(t.metadata) for t in tensors]

    with h5py.File("{}/{}.hdf5".format(conf["tensors_path"], sample), "w") as f:
        data_dset= f.create_dataset("data", data=raw_data, compression="gzip")
        labels_dset= f.create_dataset("labels", data=labels, compression="gzip")
        metadata_dset= f.create_dataset("metadata", data=metadata, compression="gzip", dtype=dt)
    logging.info("done saving {} to hdf5".format(sample))


def par_process_sample(sample):
    process_sample(global_conf, sample)
    return 0


def main():
    # The global variables are ugly, but we need them to be global so that they
    # can de passed to multi_processing
    global global_conf

    if len(sys.argv) != 3:
        print("Please provide 2 arguments: configuration.json whitelist")
        sys.exit(1)

    utils.set_logging()

    conf_fname = sys.argv[1]
    whitelist_fname = sys.argv[2]

    global_conf = utils.get_json(conf_fname)
    samples = get_whitelist(whitelist_fname)

    samples_size = len(samples)
    logging.info("There is a total of {} samples to process".format(samples_size))

    n_threads = global_conf.get("n_threads", 1)

    if n_threads > 1:
        p = Pool(n_threads)
        p.map(par_process_sample, samples)
    else:
        for sample in samples:
            process_sample(global_conf, sample)

    logging.info("Finished generating tensors")
    sys.exit(0)


if __name__ == "__main__":
    main()
