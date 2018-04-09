import numpy as np
import pysam

class RefRead(object):
    #TODO(gamazeps) make a superclass for SAM and FASTA seq like that
    def __init__(self, chrom, start, stop, seq):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.seq = seq
        assert(len(self.seq) == (self.stop - self.start))


class RefSeq(object):
    #TODO(gamazeps) make a superclass for SAM and FASTA pairs
    def __init__(self, fname):
        """
        Used for creating the reference sequence from FASTA format given a
        file.
        """
        with open(fname, "r") as f:
            lines = [l.strip() for l in f.readlines()]
            assert(len(lines) == 4)

            meta = lines[0].split(' ')
            assert(len(meta) == 4)
            seq = lines[1]
            self.first = RefRead(meta[1], int(meta[2]), int(meta[3]), seq=seq)

            meta = lines[2].split(' ')
            assert(len(meta) == 4)
            seq = lines[3]
            self.second = RefRead(meta[1], int(meta[2]), int(meta[3]), seq=seq)

    def leftmost(self):
        return self.first.pos

    def align(self, read):
        return [c for c in read.seq]

    def extract_range(self, begin, end, dummy=' '):
        assert(begin < end)
        res = np.full((end - begin), dummy)
        for read in [self.first, self.second]:
            if read:
                if read.stop < begin:
                    continue
                elif read.start > end:
                    continue
                elif read.stop == begin:
                    res[0] = self.align(read)[-1]
                elif read.start == end:
                    res[-1] = self.align(read)[0]
                elif read.start >= begin and read.stop <= end:
                    res[read.start - begin: read.stop - begin] = self.align(read)
                    assert(np.array_equal(res[read.start - begin: read.stop - begin],
                                          self.align(read)))
                elif read.start < begin and read.stop < end:
                    res[: read.stop - begin] = self.align(read)[-(read.stop - begin):]
                    assert(np.array_equal(res[:read.stop - begin],
                                          self.align(read)[-(read.stop - begin):]))
                elif read.start > begin and read.stop > end:
                    res[read.start - begin:] = self.align(read)[:-(read.stop - end)]
                    assert(np.array_equal(res[read.start - begin:],
                                          self.align(read)[:-(read.stop - end)]))
                elif read.start < begin and read.stop > end:
                    res[:] = self.align(read)[begin - read.start: -(read.stop - end)]
                    assert(np.array_equal(res[:],
                                          self.align(read)[begin - read.start: -(read.stop - end)]))
                else:
                    assert False, "unreachable"
        return res


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
        if read.is_read1:
            assert(self.first is None)
            self.first = read
        elif read.is_read2:
            assert(self.second is None)
            self.second = read
        else:
            assert False, "adding an alone read to a pair " + str(read.flag)

    def paired(self):
        return not (self.first is None or self.second is None)

    def align(self, read):
        return [c for c in read.seq]

    def extract_range(self, begin, end, dummy=' '):
        assert(begin < end)
        res = np.full((end - begin), dummy)
        for read in [self.first, self.second]:
            if read:
                # FIXME(gamazeps): we need to use the CIGAR stuff
                read_start = read.reference_start
                read_stop = read_start + len(read.seq)
                if read_stop < begin:
                    continue
                elif read_start > end:
                    continue
                elif read_stop == begin:
                    res[0] = self.align(read)[-1]
                elif read_start == end:
                    res[-1] = self.align(read)[0]
                elif read_start >= begin and read_stop <= end:
                    res[read_start - begin: read_stop - begin] = self.align(read)
                    assert(np.array_equal(res[read_start - begin: read_stop - begin],
                                          self.align(read)))
                elif read_start < begin and read_stop < end:
                    res[: read_stop - begin] = self.align(read)[-(read_stop - begin):]
                    assert(np.array_equal(res[:read_stop - begin],
                                          self.align(read)[-(read_stop - begin):]))
                elif read_start > begin and read_stop > end:
                    res[read_start - begin:] = self.align(read)[:-(read_stop - end)]
                    assert(np.array_equal(res[read_start - begin:],
                                          self.align(read)[:-(read_stop - end)]))
                elif read_start < begin and read_stop > end:
                    res[:] = self.align(read)[begin - read_start: -(read_stop - end)]
                    assert(np.array_equal(res[:],
                                          self.align(read)[begin - read_start: -(read_stop - end)]))
                else:
                    assert False, "unreachable"
        return res
