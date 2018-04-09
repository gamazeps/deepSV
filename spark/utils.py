import json
import logging
import os


def get_json(fname):
    """
    Gets the metadata from the json file
    """
    with open(fname, "r") as f:
        metadata = json.load(f)
    return metadata


def set_logging(fname=None):
    hostname = os.uname().nodename
    logging.basicConfig(level=logging.INFO,
                        filename=fname,
                        format='%(asctime)s %(levelname)-8s {} -- %(message)s'.format(hostname),
                        datefmt='%m/%d/%Y %H:%M:%S')
