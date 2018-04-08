import json
import logging


def get_json(fname):
    """
    Gets the metadata from the json file
    """
    with open(fname, "r") as f:
        metadata = json.load(f)
    return metadata


def set_logging(fname=None):
    logging.basicConfig(level=logging.INFO,
                        filename=fname,
                        format='%(asctime)s %(levelname)s -- %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')
