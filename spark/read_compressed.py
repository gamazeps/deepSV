import h5py
import time
import random
import sys
import logging

import utils

def extract_batch(fname, batchsize, iters):
    start = time.time()
    for i in range(iters):
        with h5py.File(fname, 'r') as f:
            (samples, rows, cols, channels) = f['data'].shape
            r = random.randint(0, samples - batchsize)
            _dummy = f['data'][r: r+batchsize]
    stop = time.time()
    logging.info('{} tensors, {} batchsize -- avg: {}, avg_per_tensor {}'
                 .format(samples,
                         batchsize,
                         (stop - start) / iters,
                         (stop - start) / (iters * batchsize)))
    return stop - start


def grid_search(conf):
    fnames = conf['fnames']
    batches = [16, 32, 64, 128, 256]

    for fname in fnames:
        for batch in batches:
            extract_batch(fname, batch, 3)


def main():
    utils.set_logging()

    if len(sys.argv) != 2:
        print("Please provide 1 argument: configuration.json")
        sys.exit(1)

    conf = utils.get_json(sys.argv[1])

    grid_search(conf)

if __name__ == '__main__':
    main()
