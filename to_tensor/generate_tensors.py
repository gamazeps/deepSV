import glob
from multiprocessing import Pool
import cPickle
import sys
import random

from read_pairs import SamRead, ReadPair, RefSeq
from utils import get_json
from deepsv_tensor import TensorEncoder, DeepSVTensor

global_conf = None

def find_variant_files(path, sample):
    fnames = glob.glob("{}/{}/*sam".format(path, sample))
    return [fname[:-4] for fname in fnames]


def build_tensor(basename, reads_limit=150):
    metadata = get_json(basename + ".json")
    record = metadata["record"]
    read_pairs = build_read_pairs(basename + ".sam", reads_limit)
    ref = RefSeq(basename + ".fa")
    encoder = TensorEncoder(n_channels=8, sam_channels=4, ref_channels=4)

    tensor = DeepSVTensor(encoder=encoder, metadata=record, pairs_capacity=reads_limit)
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


def process_variant(fname, index=None, draw=False):
    tensor = build_tensor(fname)
    if index:
        print(index)
    if draw:
        tensor.dummy_image(fname + ".png", draw_bp=True)
        print(index, fname + ".png")
    return tensor


def process_sample(conf, sample):
    names = find_variant_files(conf["reads_path"], sample)
    tensors = [process_variant(fname) for fname in names]

    print("start pickling")
    with open("{}/{}.pckl".format(conf["tensors_path"], sample), "wb") as f:
        cPickle.dump(tensors, f, cPickle.HIGHEST_PROTOCOL)


def par_process_sample(pair):
    process_sample(global_conf, pair[1])
    print("processed {}th sample, {}".format(pair[0], pair[1]))
    return 0


def main():
    global global_conf

    if len(sys.argv) != 3:
        print("Please provide 2 arguments: configuration.json whitelist")
        sys.exit(1)

    conf_fname = sys.argv[1]
    whitelist_fname = sys.argv[2]

    global_conf = get_json(conf_fname)
    samples = get_whitelist(whitelist_fname)

    samples_size = len(samples)
    print("There is a total of {} reads to process".format(samples_size))

    n_threads = global_conf.get("n_threads", 1)

    if n_threads > 1:
        p = Pool(n_threads)
        p.map(par_process_sample, enumerate(samples))
    else:
        for i, sample in enumerate(samples):
            print("processed {}th sample, {}".format(i, sample))
            process_sample(global_conf, sample)

    sys.exit(0)


if __name__ == "__main__":
    main()
