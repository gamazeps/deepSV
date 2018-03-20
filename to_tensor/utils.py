import json


def get_json(fname):
    """
    Gets the metadata from the json file
    """
    with open(fname, "r") as f:
        metadata = json.load(f)
    return metadata
