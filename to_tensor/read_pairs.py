import numpy as np

class RefRead(object):
    #TODO(gamazeps) make a superclass for SAM and FASTA seq like that
    def __init__(self, chrom, begin, end, seq=""):
        self.chrom = chrom
        self.pos = begin - 1 # FIXME(gamazeps): this is an ugly hack
        self.end = end
        self.seq = seq

    def align(self):
        return [c for c in self.seq]

    def append_seq(self, seq):
        self.seq += seq


class RefSeq(object):
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
                    chrom = parts[0]
                    (beg, end) = parts[1].split("-")
                    curr = RefRead(chrom, int(beg), int(end))
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


class SamRead(object):
    def __init__(self, sam_record):
        """
        TODO(gamazeps): this is awfull
        """
        tokens = sam_record.split('\t')
        self.qname = tokens[0]
        self.flag = int(tokens[1])
        self.rname = tokens[2]
        self.pos = int(tokens[3])
        self.mapq = tokens[4]
        self.cigar = tokens[5]
        self.rnext = tokens[6]
        self.pnext = tokens[7]
        self.tlen = tokens[8]
        self.seq = tokens[9]
        self.qual = tokens[10]

        self.size = len(self.seq)
        self.end = self.pos + self.size

    def align(self):
        # TODO(gamazeps): do the CIGAR stuff here.
        return [c for c in self.seq]


class ReadPair(object):
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
            assert(self.first is None)
            self.first = read
        elif read.flag & 128:
            # second red in the readpair
            assert(self.second is None)
            self.second = read
        else:
            assert False, "adding an alone read to a pair " + str(read.flag)

    def paired(self):
        return not (self.first is None or self.second is None)

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
