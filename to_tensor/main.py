import glob
from PIL import Image, ImageDraw
import cPickle
import json
import numpy as np

median_read_size = 101
median_insert_size = 500
breakpoint_window = 2 * (median_read_size + median_insert_size)
split_marker = 10
image_w = 2 * breakpoint_window + split_marker


def find_sam_files():
    return glob.glob("../data/supporting_reads/NA12878/*sam")


class RefRead:
    def __init__(self, chr, begin, end, seq=""):
        self.chr = chr
        self.begin = begin
        self.end = end
        self.seq = seq

    def append_seq(self, seq):
        self.seq += seq


class RefSeq:
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
                    curr = RefRead(chr, beg, end)
                else:
                    curr.append_seq(l.strip())
            res.append(curr)
        assert(len(res) == 2)
        self.first = res[0]
        self.second = res[1]

    def leftmost(self):
        return self.first.begin


def get_pos(pos, variant_size, split_windows):
    """
    Helper fucntion for finding the location of a base in the tensor
    """
    # We deal with case where we need to split the reads
    if split_windows:
        if pos > (variant_size - breakpoint_window):
            pos = breakpoint_window + split_marker + pos - (variant_size - breakpoint_window)
        else:
            pass
    else:
        pos += (image_w - variant_size) / 2
    return pos


def draw_sam(fname):
    f = open(fname + ".sam", "r")
    content = [sam_parser(record) for record in f]

    if len(content) is 0:
        return # we need to be careful of empty files

    origin = content[0]["POS"]
    end    = content[-1]["POS"] + len(content[-1]["SEQ"])

    reads_id = set()
    for r in content:
        reads_id.add(r["QNAME"])

    ref = read_fa(fname + ".fa")
    content = ref + content

    l = len(reads_id)

    #print(origin, end, l, len(content), len(content) - l, fname + ".png")

    img  = Image.new("RGB", (image_w, l), (0, 0, 0))
    draw = ImageDraw.Draw(img, "RGB")

    values = {
                "A": (54, 117, 177),
                "C": (241, 133, 39),
                "G": (79, 158, 57),
                "T": (199, 56, 44),
                "N": (127, 127, 127),
             }

    row = 0
    index = dict()
    variant_size = end - origin
    split_windows = variant_size > image_w
    for read in content:
        # Needed for pairing read pairs together
        j = row
        if read["QNAME"] in index:
            j = index[read["QNAME"]]
        else:
            index[read["QNAME"]] = row
            row += 1

        pos = get_pos(read["POS"] - origin, variant_size, split_windows)

        for i, base in enumerate(read["SEQ"]):
            draw.point((pos + i, j), values[base])

    # We draw the split marker
    if split_windows:
        for i in range(breakpoint_window, breakpoint_window + split_marker):
            for j in range(0, l):
                draw.point((i, j), (255, 0, 0))

    img.save(fname + ".png")


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


class DeepSVTensor:
    n_channels = 8

    def __init__(self, size, pairs_capacity, metadata):
        self.metadata = metadata
        self.size = size
        self.pairs_capacity = pairs_capacity
        self.tensor = np.zeros((size, pairs_capacity, self.n_channels))
        self.n_pairs = 0

    def insert_read_pair(self, rp):
        self.n_pairs += 1
        if self.n_pairs >= self.pairs_capacity:
            print("Too many pairs in variant: " + self.metadata["id"] )
            pass

    def insert_ref(self, ref):
        pass


def build_tensor(basename):
    metadata = get_metadata(basename + ".json")
    record = metadata["record"]
    begin = record["pos"]
    end = record["info"]["END"]["END"]
    read_pairs = build_read_pairs(basename + ".sam")
    ref = RefSeq(basename + ".fa")

    tensor = DeepSVTensor(size = 20, pairs_capacity = 200, metadata = record)
    tensor.insert_ref(ref)
    for pair in read_pairs:
        tensor.insert_read_pair(pair)

    return tensor


def build_read_pairs(fname):
    with open(fname, "r") as f:
        records = dict()
        for record in f:
            read = SamRead(record)
            if not (read.flag & int("0xc0", base=16)):
                continue
            if read.qname not in records:
                records[read.qname] = ReadPair(read.qname)
            records[read.qname].add_read(read)
    return sorted([v for (k, v) in records.iteritems()],
                  key=lambda rp: rp.leftmost())


def get_metadata(fname):
    """
    Gets the metadata from the json file
    """
    with open(fname, "r") as f:
        metadata = json.load(f)
    return metadata


if __name__ == "__main__":
    names = find_sam_files()
    for fname in names[:10]:
        fname = fname[:-len(".sam")]
        build_tensor(fname)
