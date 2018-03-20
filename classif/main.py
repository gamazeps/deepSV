import numpy as np
import sklearn as sk
import sys
import json
import cPickle


def fetch_samples(path):
    print("files are located in: {}".format(path))
    with open(path + "/NA12878.pckl", "rb") as f:
        return cPickle.load(f)


def get_json(fname):
    """
    Gets the metadata from the json file
    """
    with open(fname, "r") as f:
        metadata = json.load(f)
    return metadata


def main():
    if len(sys.argv) != 2:
        print("Please provide 1 arguments: configuration.json")
        sys.exit(1)

    conf_fname = sys.argv[1]
    conf = get_json(conf_fname)

    data = fetch_samples(conf["tensors_path"])

    sys.exit(0)


if __name__ == "__main__":
    main()
