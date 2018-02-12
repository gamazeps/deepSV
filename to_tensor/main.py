import glob
from PIL import Image, ImageDraw

median_read_size = 101
median_insert_size = 400
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


def draw_sam(fname):
    f = open(fname, "r")
    content = [ugly_parse(record) for record in f]

    if len(content) is 0:
        return # we need to be careful of empty files

    reads_id = set()
    for r in content:
        reads_id.add(r["QNAME"])

    origin = content[0]["POS"]
    end    = content[-1]["POS"] + len(content[-1]["SEQ"])

    l = len(reads_id)

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
        j = row
        pos = read["POS"] - origin

        # We deal with case where we need to split the reads
        if split_windows:
            if pos > (variant_size - breakpoint_window):
                pos = breakpoint_window + split_marker + pos - (variant_size - breakpoint_window)
                #assert(pos + len(read["SEQ"]) < variant_size)
            else:
                #assert(pos + len(read["SEQ"]) < breakpoint_window)
                pass
        else:
            pos += (image_w - variant_size) / 2

        # Needed for pairing read pairs together
        if read["QNAME"] in index:
            j = index[read["QNAME"]]
        else:
            index[read["QNAME"]] = row
            row += 1

        for i, base in enumerate(read["SEQ"]):
            draw.point((pos + i, j), values[base])

    # We draw the split marker
    if split_windows:
        for i in range(breakpoint_window, breakpoint_window + split_marker):
            for j in range(0, l):
                draw.point((i, j), (255, 0, 0))


    img.save(fname + ".png")


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
