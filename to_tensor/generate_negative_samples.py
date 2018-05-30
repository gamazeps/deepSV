import argparse
import collections
import random
import json

import numpy as np
import pysam

# Taken from wikipedia
chr_len = {
    "1" :  248956422,
    "2" :  242193529,
    "3" :  198295559,
    "4" :  190214555,
    "5" :  181538259,
    "6" :  170805979,
    "7" :  159345973,
    "8" :  145138636,
    "9" :  138394717,
    "10":  133797422,
    "11":  135086622,
    "12":  133275309,
    "13":  114364328,
    "14":  107043718,
    "15":  101991189,
    "16":  90338345 ,
    "17":  83257441 ,
    "18":  80373285 ,
    "19":  58617616 ,
    "20":  64444167 ,
    "21":  46709983 ,
    "22":  50818468 ,
    "X" :  156040895,
    "Y" :  57227415 ,
}

class NegativeRecord():
    def __init__(self, id, chrom, start, stop, info):
        self.id = id
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.info = info

def keep_record(record):
    # We want to remove variants with uncertain locations
    return not ("CIEND" in record.info.keys() or "CIPOS" in record.info.keys())

def main():
    parser = argparse.ArgumentParser(
        description='Tool for genrating negative samples')
    parser.add_argument("--vcf-path", type=str, required=True)
    parser.add_argument("--negative-ratio", type=float, required=True)
    parser.add_argument("--from-pos-region", type=float, required=True)
    parser.add_argument("--outpath", type=str, required=True)
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf_path, 'rb')
    records = [record for record in vcf.fetch() if keep_record(record)]


    # statistics for the variants
    individuals_for_region = collections.defaultdict(set)
    region_count = collections.defaultdict(int)
    variants_per_individual = collections.defaultdict(int)
    region_info = dict()

    for record in records:
        individuals_for_region[record.info["REGIONID"][0]].add(record.info["SAMPLE"])
        region_count[record.info["REGIONID"][0]] += 1
        variants_per_individual[record.info["SAMPLE"]]+=1
        region_info[record.info["REGIONID"][0]] = record

    random.seed(0) # we fix the seed for reproducibility

    res = list()

    for region in region_count.keys():
        local_individuals = set(variants_per_individual)
        local_individuals -= individuals_for_region[region]
        n_samples = min(
            int(region_count[region] * args.negative_ratio * args.from_pos_region),
            len(local_individuals))

        samples = random.sample(list(local_individuals), n_samples)
        record = region_info[region]

        for i, sample in enumerate(samples):
            neg_record = {
                "id": "fake-{}-{}".format(region, i),
                "chrom": record.chrom,
                "start": record.start,
                "stop": record.stop,
                "info": {"SAMPLE": sample, "SVTYPE": "FAKE"}
            }
            res.append(neg_record)

    with open(args.outpath, "w") as f:
        f.write(json.dumps(res))

if __name__ == "__main__":
    main()
