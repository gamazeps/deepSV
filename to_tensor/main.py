import glob
from multiprocessing import Pool
from PIL import Image, ImageDraw
import cPickle
import json
import numpy as np
import sys
import random

median_read_size = 101
median_insert_size = 400
breakpoint_window = 2 * (median_read_size + median_insert_size)
split_marker = 10
full_window = 2 * breakpoint_window + split_marker


def find_variant_files(path, sample):
    fnames = glob.glob("{}/{}/*sam".format(path, sample))
    return [fname[:-4] for fname in fnames]


class RefRead:
    #TODO(gamazeps) make a superclass for SAM and FASTA seq like that
    def __init__(self, chr, begin, end, seq=""):
        self.chr = chr
        self.pos = begin - 1 # FIXME(gamazeps): this is an ugly hack
        self.end = end
        self.seq = seq

    def align(self):
        return [c for c in self.seq]

    def append_seq(self, seq):
        self.seq += seq


class RefSeq:
    #TODO(gamazeps) make a superclass for SAM and FASTA pairs
    def __init__(self, fname):
        """
        Used for creating the reference sequence from FASTA format given a
        file.
        """
        with open(fname, "r") as f:
            res = []
            curr = None
            for l in f:
                if l[0] == ">":
                    if curr is not None:
                        res.append(curr)
                    parts = l[1:].strip().split(":")
                    chr = parts[0]
                    (beg, end) = parts[1].split("-")
                    curr = RefRead(chr, int(beg), int(end))
                else:
                    curr.append_seq(l.strip())
            assert(len(curr.seq) == (curr.end - curr.pos))
            res.append(curr)
        assert(len(res) == 2)
        self.first = res[0]
        self.second = res[1]

    def leftmost(self):
        return self.first.pos

    def extract_range(self, begin, end, dummy=' '):
        assert(begin < end)
        res = np.full((end - begin), dummy)
        for read in [self.first, self.second]:
            if read:
                if read.end < begin:
                    continue
                elif read.pos > end:
                    continue
                elif read.end == begin:
                    res[0] = read.align()[-1]
                elif read.pos == end:
                    res[-1] = read.align()[0]
                elif read.pos >= begin and read.end <= end:
                    res[read.pos - begin: read.end - begin] = read.align()
                    assert(np.array_equal(res[read.pos - begin: read.end - begin],
                                          read.align()))
                elif read.pos < begin and read.end < end:
                    res[: read.end - begin] = read.align()[-(read.end - begin):]
                    assert(np.array_equal(res[:read.end - begin],
                        read.align()[-(read.end - begin):]))
                elif read.pos > begin and read.end > end:
                    res[read.pos - begin:] = read.align()[:-(read.end - end)]
                    assert(np.array_equal(res[read.pos - begin:],
                                          read.align()[:-(read.end - end)]))
                elif read.pos < begin and read.end > end:
                    res[:] = read.align()[begin - read.pos: -(read.end - end)]
                    assert(np.array_equal(res[:],
                                          read.align()[begin - read.pos: -(read.end - end)]))
                else:
                    assert False, "unreachable"
        return res


class SamRead:
    def __init__(self, sam_record):
        """
        TODO(gamazeps): this is awfull
        """
        tokens = sam_record.split('\t')
        self.qname = tokens[0]
        self.flag  = int(tokens[1])
        self.rname = tokens[2]
        self.pos   = int(tokens[3])
        self.mapq  = tokens[4]
        self.cigar = tokens[5]
        self.rnext = tokens[6]
        self.pnext = tokens[7]
        self.tlen  = tokens[8]
        self.seq   = tokens[9]
        self.qual  = tokens[10]

        self.size = len(self.seq)
        self.end = self.pos + self.size

    def align(self):
        # TODO(gamazeps): do the CIGAR stuff here.
        return [c for c in self.seq]


class TensorEncoder:
    def __init__(self, n_channels=4, sam_channels=4, ref_channels=0):
        assert(n_channels == (sam_channels + ref_channels))
        self.n_channels = n_channels
        self.sam_channels = sam_channels
        self.ref_channels = ref_channels

    def encode_sam(self, read):
        """
        hardcoded encoding, this is bad
        """
        encoding = {
                " ": [0, 0, 0, 0],
                "N": [0, 0, 0, 0],
                "A": [1, 0, 0, 0],
                "C": [0, 1, 0, 0],
                "G": [0, 0, 1, 0],
                "T": [0, 0, 0, 1],
        }

        (x,) = read.shape
        data = np.full((x, self.sam_channels), 0)
        for i in range(0, x):
            data[i] = encoding[read[i]]
        return data

    def encode_ref(self, read):
        """
        hardcoded encoding, this is bad
        NOTE: copy of `encode_sam` for now, but sam will soon have more info (qual for ex)
        """
        encoding = {
                " ": [0, 0, 0, 0],
                "N": [0, 0, 0, 0],
                "A": [1, 0, 0, 0],
                "C": [0, 1, 0, 0],
                "G": [0, 0, 1, 0],
                "T": [0, 0, 0, 1],
        }

        (x,) = read.shape
        data = np.full((x, self.ref_channels), 0)
        for i in range(0, x):
            data[i] = encoding[read[i]]
        return data

    def decode_sam(self, encoding, dummy=' '):
        """
        hardcoded encoding, this is bad
        """
        (x, y) = encoding.shape
        assert(y == self.n_channels)
        res = np.full(x, dummy)
        for i in range(0, x):
            # TODO(gamazeps): this is probably horribly slow
            if encoding[i][0] == 1:
                res[i] = "A"
            elif encoding[i][1] == 1:
                res[i] = "C"
            elif encoding[i][2] == 1:
                res[i] = "G"
            elif encoding[i][3] == 1:
                res[i] = "T"
        return res


class ReadPair:
    def __init__(self, name):
        self.name = name
        self.first = None
        self.second = None

    def leftmost(self):
        if self.first:
            return self.first.pos
        elif self.second:
            return self.second.pos
        else:
            assert(False)

    def add_read(self, read):
        if read.flag & 64:
            # first read of the pair
            assert(self.first == None)
            self.first = read
        elif read.flag & 128:
            # second red in the readpair
            assert(self.second == None)
            self.second = read
        else:
            assert False, "adding an alone read to a pair " + str(read.flag)

    def paired(self):
        return not (self.first == None or self.second == None)

    def extract_range(self, begin, end, dummy=' '):
        assert(begin < end)
        res = np.full((end - begin), dummy)
        for read in [self.first, self.second]:
            if read:
                if read.end < begin:
                    continue
                elif read.pos > end:
                    continue
                elif read.end == begin:
                    res[0] = read.align()[-1]
                elif read.pos == end:
                    res[-1] = read.align()[0]
                elif read.pos >= begin and read.end <= end:
                    res[read.pos - begin: read.end - begin] = read.align()
                    assert(np.array_equal(res[read.pos - begin: read.end - begin],
                                          read.align()))
                elif read.pos < begin and read.end < end:
                    res[: read.end - begin] = read.align()[-(read.end - begin):]
                    assert(np.array_equal(res[:read.end - begin],
                        read.align()[-(read.end - begin):]))
                elif read.pos > begin and read.end > end:
                    res[read.pos - begin:] = read.align()[:-(read.end - end)]
                    assert(np.array_equal(res[read.pos - begin:],
                                          read.align()[:-(read.end - end)]))
                elif read.pos < begin and read.end > end:
                    res[:] = read.align()[begin - read.pos: -(read.end - end)]
                    assert(np.array_equal(res[:],
                                          read.align()[begin - read.pos: -(read.end - end)]))
                else:
                    assert False, "unreachable"
        return res


class DeepSVTensor:
    dummy = ' '

    def __init__(self, encoder, metadata, pairs_capacity):
        self.encoder = encoder
        self.metadata = metadata
        self.begin = self.metadata["pos"] - (breakpoint_window // 2)
        self.end = self.metadata["info"]["END"]["END"]
        self.size = self.end - self.begin
        self.pairs_capacity = pairs_capacity
        self.tensor = np.full((pairs_capacity, full_window, encoder.n_channels), 0)
        self.n_pairs = 0
        self.is_split = self.size > breakpoint_window
        self.label = self.metadata["alt"][0]

    def insert_read_pair(self, rp):
        if self.n_pairs >= self.pairs_capacity:
            print("Too many pairs in variant: " + self.metadata["id"] )
            self.n_pairs += 1
            return
        if self.is_split:
            self.tensor[self.n_pairs, :breakpoint_window, :self.encoder.sam_channels] = self.encoder.encode_sam(
                    rp.extract_range(
                        self.begin - (breakpoint_window // 2),
                        self.begin + (breakpoint_window // 2),
                        dummy=self.dummy
                    )
            )
            self.tensor[self.n_pairs, breakpoint_window + split_marker:, :self.encoder.sam_channels] = self.encoder.encode_sam(
                    rp.extract_range(
                        self.end - (breakpoint_window // 2),
                        self.end + (breakpoint_window // 2),
                        dummy=self.dummy
                    )
            )
        else:
            l_center_pad = (full_window - (self.size + breakpoint_window) + 1) // 2
            r_center_pad = (full_window - (self.size + breakpoint_window)) // 2
            self.tensor[self.n_pairs, l_center_pad: -r_center_pad, :self.encoder.sam_channels] = self.encoder.encode_sam(
                    rp.extract_range(
                        self.begin - (breakpoint_window // 2),
                        self.end + (breakpoint_window // 2),
                        dummy=self.dummy
                    )
            )

        self.n_pairs += 1

    def insert_ref(self, ref):
        if self.is_split:
            self.tensor[:, :breakpoint_window, self.encoder.sam_channels:
                    self.encoder.sam_channels + self.encoder.ref_channels] = self.encoder.encode_ref(
                    ref.extract_range(
                        self.begin - (breakpoint_window // 2),
                        self.begin + (breakpoint_window // 2),
                        dummy=self.dummy
                    )
            )
            self.tensor[:, breakpoint_window + split_marker:, self.encoder.sam_channels:
                    self.encoder.sam_channels + self.encoder.ref_channels] = self.encoder.encode_ref(
                    ref.extract_range(
                        self.end - (breakpoint_window // 2),
                        self.end + (breakpoint_window // 2),
                        dummy=self.dummy
                    )
            )
        else:
            l_center_pad = (full_window - (self.size + breakpoint_window) + 1) // 2
            r_center_pad = (full_window - (self.size + breakpoint_window)) // 2
            self.tensor[:, l_center_pad: -r_center_pad, self.encoder.sam_channels:
                    self.encoder.sam_channels + self.encoder.ref_channels] = self.encoder.encode_ref(
                    ref.extract_range(
                        self.begin - (breakpoint_window // 2),
                        self.end + (breakpoint_window // 2),
                        dummy=self.dummy
                    )
            )


    def dummy_image(self, fname, draw_bp=False):
        if self.pairs_capacity == 0:
            return

        img  = Image.new("RGB", (full_window, self.pairs_capacity), (0, 0, 0))
        draw = ImageDraw.Draw(img, "RGB")
        values = {
                    "A": (54,  117, 177),
                    "C": (241, 133, 39),
                    "G": (79,  158, 57),
                    "T": (199, 56,  44),
                    "N": (127, 127, 127),
                    self.dummy: (0,   0,   0)
                 }

        for i in range(0, self.pairs_capacity):
            read = self.encoder.decode_sam(self.tensor[i], dummy=self.dummy)
            for j in range(0, full_window):
                draw.point((j, i), values[read[j]])

        if self.is_split:
            for i in range(0, self.pairs_capacity):
                for j in range(breakpoint_window, breakpoint_window + split_marker):
                    draw.point((j, i), (255, 0, 0))

        if draw_bp:
            for i in range(0, self.pairs_capacity):
                draw.point((breakpoint_window / 2, i), (255, 255, 255))
                draw.point((full_window - (breakpoint_window / 2), i), (255, 255, 255))

        img.save(fname)


def build_tensor(basename, reads_limit=150):
    metadata = get_json(basename + ".json")
    record = metadata["record"]
    read_pairs = build_read_pairs(basename + ".sam", reads_limit)
    ref = RefSeq(basename + ".fa")
    encoder = TensorEncoder(n_channels=8, sam_channels=4, ref_channels=4)

    tensor = DeepSVTensor(encoder=encoder, metadata=record, pairs_capacity=reads_limit)
    tensor.insert_ref(ref)
    for pair in read_pairs:
        tensor.insert_read_pair(pair)

    return tensor


def build_read_pairs(fname, sample_limit):
    with open(fname, "r") as f:
        records = dict()
        for record in f:
            read = SamRead(record)
            if not (read.flag & int("0xc0", base=16)):
                continue
            if read.qname not in records:
                records[read.qname] = ReadPair(read.qname)
            records[read.qname].add_read(read)

        # Here we sample up to `sample_limit` reads at random
        records_list = [v for (k, v) in records.iteritems()]
        n_records = len(records_list)
        sampled_indices = sorted(random.sample(range(n_records), min(n_records, sample_limit)))
        sampled_records_list = [records_list[i] for i in sampled_indices]

        return sorted(sampled_records_list, key=lambda rp: rp.leftmost())


def get_json(fname):
    """
    Gets the metadata from the json file
    """
    with open(fname, "r") as f:
        metadata = json.load(f)
    return metadata


def get_samples(fname):
    with open(fname, "r") as f:
        return [l.strip() for l in f]


def process_variant(fname, index=None, draw=False):
    tensor = build_tensor(fname)
    if index:
        print(index)
    if draw:
        tensor.dummy_image(fname + ".png", draw_bp=True)
        print(index, fname + ".png")
    return tensor


def process_sample(conf, sample):
    names = find_variant_files(conf["reads_path"], sample)
    for i, fname in enumerate(names[:]):
        _ = process_variant(fname, index=i, draw=False)


def main():
    if len(sys.argv) != 3:
        print("Please provide 2 arguments: configuration.json whitelist")
        sys.exit(1)
    conf_fname = sys.argv[1]
    whitelist_fname = sys.argv[2]

    conf = get_json(conf_fname)
    samples = get_samples(whitelist_fname)

    for sample in samples:
        process_sample(conf, sample)

    sys.exit(0)


if __name__ == "__main__":
    main()
