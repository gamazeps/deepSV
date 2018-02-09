import glob
from PIL import Image, ImageDraw

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

    origin = content[0]["POS"]
    end    = content[-1]["POS"] + len(content[-1]["SEQ"])
    w      = end - origin
    l      = len(content)

    img  = Image.new("RGB", (w, l), (0, 0, 0))
    draw = ImageDraw.Draw(img, "RGB")

    values = {
                "A": (54, 117, 177),
                "C": (241, 133, 39),
                "G": (79, 158, 57),
                "T": (199, 56, 44),
             }

    for j, read in enumerate(content):
        pos = read["POS"] - origin
        for i, base in enumerate(read["SEQ"]):
            draw.point((pos + i, j), values[base])

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
