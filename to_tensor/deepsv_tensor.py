"""
Module for defining the DeepSVTensor class.
"""

import numpy as np
from PIL import Image, ImageDraw

MEDIAN_READ_SIZE = 101
MEDIAN_INSERT_SIZE = 400
BREAKPOINT_WINDOW = 2 * (MEDIAN_READ_SIZE + MEDIAN_INSERT_SIZE)
SPLIT_MARKER = 10
FULL_WINDOW = 2 * BREAKPOINT_WINDOW + SPLIT_MARKER

class TensorEncoder(object):
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


class DeepSVTensor(object):
    # pylint: disable=too-many-instance-attributes
    dummy = ' '

    def __init__(self, encoder, metadata, pairs_capacity):
        self.encoder = encoder
        self.metadata = metadata
        self.begin = self.metadata["start"] - (BREAKPOINT_WINDOW // 2)
        self.end = self.metadata["stop"]
        self.size = self.end - self.begin
        self.pairs_capacity = pairs_capacity
        self.tensor = np.full((pairs_capacity, FULL_WINDOW, encoder.n_channels), 0)
        self.n_pairs = 0
        self.is_split = self.size > BREAKPOINT_WINDOW

    def insert_read_pair(self, rp):
        """
        Adds a read pair to the tensor.
        """
        if self.n_pairs >= self.pairs_capacity:
            print("Too many pairs in variant: " + self.metadata["id"])
            self.n_pairs += 1
            return
        if self.is_split:
            self.tensor[self.n_pairs, :BREAKPOINT_WINDOW, :self.encoder.sam_channels] = (
                self.encoder.encode_sam(
                    rp.extract_range(
                        self.begin - (BREAKPOINT_WINDOW // 2),
                        self.begin + (BREAKPOINT_WINDOW // 2),
                        dummy=self.dummy
                    )
                )
            )
            self.tensor[self.n_pairs, BREAKPOINT_WINDOW + SPLIT_MARKER:,
                        :self.encoder.sam_channels] = (
                            self.encoder.encode_sam(
                                rp.extract_range(
                                    self.end - (BREAKPOINT_WINDOW // 2),
                                    self.end + (BREAKPOINT_WINDOW // 2),
                                    dummy=self.dummy
                                )
                            )
                        )
        else:
            l_center_pad = (FULL_WINDOW - (self.size + BREAKPOINT_WINDOW) + 1) // 2
            r_center_pad = (FULL_WINDOW - (self.size + BREAKPOINT_WINDOW)) // 2
            self.tensor[self.n_pairs, l_center_pad: -r_center_pad, :self.encoder.sam_channels] = (
                self.encoder.encode_sam(
                    rp.extract_range(
                        self.begin - (BREAKPOINT_WINDOW // 2),
                        self.end + (BREAKPOINT_WINDOW // 2),
                        dummy=self.dummy
                    )
                )
            )

        self.n_pairs += 1

    def insert_ref(self, ref):
        """
        Adds the reference sequence to the tensor.
        """
        if self.is_split:
            self.tensor[:, :BREAKPOINT_WINDOW, self.encoder.sam_channels:
                        self.encoder.sam_channels + self.encoder.ref_channels] = (
                            self.encoder.encode_ref(
                                ref.extract_range(
                                    self.begin - (BREAKPOINT_WINDOW // 2),
                                    self.begin + (BREAKPOINT_WINDOW // 2),
                                    dummy=self.dummy
                                )
                            )
                        )
            self.tensor[:, BREAKPOINT_WINDOW + SPLIT_MARKER:, self.encoder.sam_channels:
                        self.encoder.sam_channels + self.encoder.ref_channels] = (
                            self.encoder.encode_ref(
                                ref.extract_range(
                                    self.end - (BREAKPOINT_WINDOW // 2),
                                    self.end + (BREAKPOINT_WINDOW // 2),
                                    dummy=self.dummy
                                )
                            )
                        )
        else:
            l_center_pad = (FULL_WINDOW - (self.size + BREAKPOINT_WINDOW) + 1) // 2
            r_center_pad = (FULL_WINDOW - (self.size + BREAKPOINT_WINDOW)) // 2
            self.tensor[:, l_center_pad: -r_center_pad, self.encoder.sam_channels:
                        self.encoder.sam_channels + self.encoder.ref_channels] = (
                            self.encoder.encode_ref(
                                ref.extract_range(
                                    self.begin - (BREAKPOINT_WINDOW // 2),
                                    self.end + (BREAKPOINT_WINDOW // 2),
                                    dummy=self.dummy
                                )
                            )
                        )


    def dummy_image(self, fname, draw_bp=False):
        if self.pairs_capacity == 0:
            return

        img = Image.new("RGB", (FULL_WINDOW, self.pairs_capacity), (0, 0, 0))
        draw = ImageDraw.Draw(img, "RGB")
        values = {
            "A": (54, 117, 177),
            "C": (241, 133, 39),
            "G": (79, 158, 57),
            "T": (199, 56, 44),
            "N": (127, 127, 127),
            self.dummy: (0, 0, 0)
        }

        for i in range(0, self.pairs_capacity):
            read = self.encoder.decode_sam(self.tensor[i], dummy=self.dummy)
            for j in range(0, FULL_WINDOW):
                draw.point((j, i), values[read[j]])

        if self.is_split:
            for i in range(0, self.pairs_capacity):
                for j in range(BREAKPOINT_WINDOW, BREAKPOINT_WINDOW + SPLIT_MARKER):
                    draw.point((j, i), (255, 0, 0))

        if draw_bp:
            for i in range(0, self.pairs_capacity):
                draw.point((BREAKPOINT_WINDOW / 2, i), (255, 255, 255))
                draw.point((FULL_WINDOW - (BREAKPOINT_WINDOW / 2), i), (255, 255, 255))

        img.save(fname)

    def label(self):
        label_to_int = {
            "DEL": 1,
            "DUP": 2,
            "INV": 3,
            "INS": 4,
            "CNV": 5
        }
        return label_to_int[self.metadata["info"]["SVTYPE"]]
