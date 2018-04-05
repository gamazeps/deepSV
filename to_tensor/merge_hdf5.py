import h5py
import glob
import logging
import sys
from multiprocessing import Pool
import os

import utils

def pair_files(fnames):
    """
    Helper function to have the samples 2 by 2
    """
    res = []
    for i in range(0, len(fnames) // 2):
        res.append([fnames[2*i], fnames[2*i + 1]])

    if len(fnames) % 2 == 1:
        res.append([fnames[-1]])

    return res


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
                                          dtype="int64", compression="lzf")
        labels_dset = root_f.create_dataset("labels", shape=(n_var,),
                                            dtype="i8", compression="lzf")
        metadata_dset = root_f.create_dataset("metadata", shape=(n_var,),
                                              dtype=dt, compression="lzf")
        curr = 0
        for fname in fnames:
            with h5py.File(fname, "r") as f:
                (curr_var,) = f["labels"].shape
                data_dset[curr: curr + curr_var] = f["data"]
                metadata_dset[curr: curr + curr_var] = f["metadata"]
                labels_dset[curr: curr + curr_var] = f["labels"]
                curr += curr_var


def par_merge(pair):
    (i, fnames) = pair
    new_name = "{}/{}.hdf5".format(out_dir, i)
    merge_samples(fnames, new_name)
    logging.info("Done merging {}".format(fnames))
    return new_name


def main():
    global out_dir

    if len(sys.argv) != 2:
        print("Please provide 1 arguments: conf.json")
        sys.exit(1)

    conf = utils.get_json(sys.argv[1])

    in_dir = conf["in_dir"]
    root_dir = conf["root_dir"]

    step = 0

    utils.set_logging()

    fnames = glob.glob("{}/*.hdf5".format(in_dir))

    if n_threads > 1:
        p = Pool(n_threads)

    while len(fnames) > 1:
        logging.info("Starting merge phase {}".format(step))
        out_dir = "{}/iter_{}".format(root_dir, step)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        paired_fnames = pair_files(fnames)

        logging.info("Will merge {} files".format(len(fnames)))

        n_threads = conf.get("n_threads", 1)

        if n_threads > 1:
            fnames = p.map(par_merge, enumerate(paired_fnames))
        else:
            fnames = [par_merge(pair) for pair in enumerate(paired_fnames)]

        logging.info("Finished merge phase {}".format(step))
        step += 1

    logging.info("Finished merging, final folder is: {}".format(out_dir))


if __name__ == "__main__":
    main()
