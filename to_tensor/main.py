import glob
from PIL import Image, ImageDraw
import cPickle

median_read_size = 101
median_insert_size = 500
breakpoint_window = 2 * (median_read_size + median_insert_size)
split_marker = 10
image_w = 2 * breakpoint_window + split_marker

def find_sam_files():
    return glob.glob("../data/supporting_reads/NA12878/*sam")

def read_sam(fname):
    f = open(fname, "r")
    content = f.readlines()

    for c in content:
        inf = ugly_parse(c)
        print(inf["POS"], len(inf["SEQ"]))

def read_fa(fname):
    f = open(fname, "r")

    res = []
    curr = None
    for l in f:
        if l[0] == ">":
            if curr is not None:
                res.append(curr)
            curr = {}
            parts = l[1:].strip().split(":")
            print("ref:" + parts[1])
            regions = parts[1].split("-")
            curr["QNAME"] = "ref"
            curr["POS"] = int(regions[0]) # we don't care about chr nor end
            curr["SEQ"] = str("")
        else:
            curr["SEQ"] += l.strip()

    res.append(curr)
    return res


def get_pos(pos, variant_size, split_windows):
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
    f = open(fname, "r")
    content = [ugly_parse(record) for record in f]

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

    print(origin, end, l, len(content), len(content) - l, fname + ".png")

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


def serialize_sam(fname):
    f = open(fname, "r")
    content = [ugly_parse(record) for record in f]

    if len(content) is 0:
        return # we need to be careful of empty files

    cPickle.dump(content, open(fname + ".pckl", "w"))
    copy = cPickle.load(open(fname + ".pckl", "r"))

    assert(copy == content)


def ugly_parse(sam):
    tokens = sam.split('\t')
    parsed = dict()

    parsed["QNAME"] = tokens[0]
    parsed["FLAG"]  = tokens[1]
    parsed["RNAME"] = tokens[2]
    parsed["POS"]   = int(tokens[3])
    parsed["MAPQ"]  = tokens[4]
    parsed["CIGAR"] = tokens[5]
    parsed["RNEXT"] = tokens[6]
    parsed["PNEXT"] = tokens[7]
    parsed["TLEN"]  = tokens[8]
    parsed["SEQ"]   = tokens[9]
    parsed["QUAL"]  = tokens[10]

    return parsed

if __name__ == "__main__":
    names = find_sam_files()
    for fname in names:
        draw_sam(fname)
        serialize_sam(fname)
