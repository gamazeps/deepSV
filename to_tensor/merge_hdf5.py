import h5py
import glob
import logging
import sys
from multiprocessing import Pool

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


def merge_samples(samples, ofile):
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
    (fnames, i) = pair
    merge_samples(fnames, "{}/{}.hdf5".format(out_dir, i))
    logging.info("Done merging {}".format(fnames))
    return None


def main():
    global out_dir

    if len(sys.argv) != 4:
        print("Please provide 3 arguments: conf.json in_dir out_dir")
        sys.exit(1)

    in_dir = sys.argv[2]
    out_dir = sys.argv[3]

    utils.set_logging()

    logging.info("Will merge {} files".format(len(fnames)))

    fnames = glob.glob("{}/*.hdf5".format(in_dir))
    paired_fnames = pair_files(fnames)

    n_threads = conf.get("n_threads", 1)

    if n_threads > 1:
        p = Pool(n_threads)
        p.map(par_merge, enumerate(paired_fnames))
    else:
        [par_merge(pair) for pair in enumerate(paired_fnames)]

    logging.info("Finished merging")


if __name__ == "__main__":
    main()
