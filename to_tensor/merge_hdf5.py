import h5py
import glob
import logging
import sys

import utils

def main():
    if len(sys.argv) != 2:
        print("Please provide 1 arguments configuration.json")
        sys.exit(1)

    utils.set_logging()

    conf = utils.get_json(sys.argv[1])
    fnames = glob.glob("{}/*.hdf5".format(conf["samples_path"]))

    logging.info("Will merge {} files".format(len(fnames)))

    # We need the total amount of variants before merging
    n_var = 0
    (window_size, reads, channels)  = (None, None, None)
    for fname in fnames:
        with h5py.File(fname, "r") as f:
            (curr_var, curr_window_size, curr_reads, curr_channels) = f["data"].shape
            logging.info("{} variants in {}".format(curr_var, fname))
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

    with h5py.File(conf["merged_tensor_path"], "w") as root_f:
        data_dset = root_f.create_dataset("data", shape=(n_var, window_size, reads, channels),
                                          dtype="int64", compression="gzip")
        labels_dset = root_f.create_dataset("labels", shape=(n_var,),
                                            dtype="int64", compression="gzip")
        metadata_dset = root_f.create_dataset("metadata", shape=(n_var,),
                                              dtype=dt, compression="gzip")
        curr = 0
        for fname in fnames:
            with h5py.File(fname, "r") as f:
                (curr_var,) = f["labels"].shape
                data_dset[curr: curr + curr_var] = f["data"]
                metadata_dset[curr: curr + curr_var] = f["metadata"]
                labels_dset[curr: curr + curr_var] = f["labels"]
                curr += curr_var
                logging.info("Done copying {} variants out of {}".format(curr, n_var))

    logging.info("Finished merging")



if __name__ == "__main__":
    main()
