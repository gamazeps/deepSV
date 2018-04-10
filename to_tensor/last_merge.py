import h5py
import glob
import logging
import sys
from multiprocessing import Pool
import os

import utils

def merge_samples(fnames, ofile):
    # We need the total amount of variants before merging
    n_var = 0
    (window_size, reads, channels)  = (None, None, None)
    for fname in fnames:
        with h5py.File(fname, "r") as f:
            (curr_var, curr_window_size, curr_reads, curr_channels) = f["data"].shape
            n_var += curr_var

            # Check for proper formatting
            if window_size is not None:
                assert(window_size == curr_window_size)
            if channels is not None:
                assert(channels == curr_channels)
            if window_size is not None:
                assert(window_size == curr_window_size)
            (window_size, reads, channels)  = (curr_window_size, curr_reads, curr_channels)

    # Needed for encoding the json metadata
    dt = h5py.special_dtype(vlen=bytes)

    with h5py.File(ofile, "w") as root_f:
        data_dset = root_f.create_dataset("data", shape=(n_var, window_size, reads, channels),
                                          dtype="u1", compression="lzf")
        labels_dset = root_f.create_dataset("labels", shape=(n_var,),
                                            dtype="u1", compression="lzf")
        metadata_dset = root_f.create_dataset("metadata", shape=(n_var,),
                                              dtype=dt, compression="lzf")
        curr = 0
        for fname in fnames:
            logging.info("merging {}".format(fname))
            with h5py.File(fname, "r") as f:
                (curr_var,) = f["labels"].shape
                step = 20000
                for i in range(0, curr_var, step)
                    local_end = min(curr_var, i + step)
                    file_end = min(curr + curr_var, curr + i + step)
                    data_dset[curr + i: file_end] = f["data"][i: local_end]
                    metadata_dset[curr + i: file_end] = f["metadata"][i: local_end]
                    labels_dset[curr + i: file_end] = f["labels"][i: local_end]
                    logging.info("Logged {} out of {} varians {}%"
                                 .format(file_end, n_var, float(file_end) / n_var))
                curr += curr_var


def main():
    if len(sys.argv) != 3:
        print("Please provide 2 arguments: conf.json log")
        sys.exit(1)

    conf = utils.get_json(sys.argv[1])

    in_dir = conf["in_dir"]
    out_dir = conf["out_dir"]
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    utils.set_logging(sys.argv[2])

    fnames = glob.glob(os.join(in_dir, "*.h5"))

    logging.info("Starting merging".format(step))
    merge_samples(fnames, os.join(out_dir, 'variants.h5'))
    logging.info("Finished merging".format(step))

if __name__ == "__main__":
    main()

