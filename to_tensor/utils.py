import json
import logging


def get_json(fname):
    """
    Gets the metadata from the json file
    """
    with open(fname, "r") as f:
        metadata = json.load(f)
    return metadata


def set_logging():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s -- %(message)s',
                        datefmt='%m/%d/%Y %H:%M:%S')
