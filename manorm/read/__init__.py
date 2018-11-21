"""Module for reads and related operation."""
import os
from collections import defaultdict
from manorm.logger import logger
from bisect import bisect_left
READ_FORMATS = ['bed', 'sam', 'bam']
class ReadFormatError(Exception):
    pass


class Reads(object):
    """Class for reads generated from next-generation sequencing."""

    def __init__(self, path=None, name=None, shift=0):
        self.path = path
        self.name = name
        if self.path and not self.name:
            self.name = os.path.splitext(os.path.basename(path))[0]
        self.pos = defaultdict(list)
        self._size = None
        if self.path:
            self.read_file(path, shift)
            self.sort()


    def __len__(self):
        return self.size

    @property
    def size(self):
        if not self._size:
            num = 0
            for chrom in self.pos:
                num += len(self.pos[chrom])
            self._size = num
        return self._size

    def sort(self):
        for chrom in self.pos:
            self.pos[chrom].sort()

    def _count_reads(self, reads, extend):
        """Count reads by binary search."""
        try:
            head = bisect_left(reads.pos[self.chrom], self.summit - extend)
            tail = bisect_left(reads.pos[self.chrom], self.summit + extend)
            return tail - head
        except KeyError:
            return 0


def load_reads(path, paired=False, shift=100):
    """Load positions of sequencing reads from input file."""
    with open(path, "r") as fin:
        for line in fin:
            fields = line.strip().split("\t")
            chrom, start, end, strand = fields[0], int(fields[1]), int(fields[2]), fields[5]
            if strand == '+':
                pos = start + shift
            elif strand == '-':
                pos = end - shift
            else:
                logger.warning("Invalid Format of Read: {}".format(line.rstrip()))
            pos[chrom].append(pos)
#
# def infer_peak_format(path):
#     """Automatically infer the format of input peak file."""
#     for format in PEAK_FORMATS:
#         logger.debug("Matching with {!r} format".format(format))
#         parser = peak_parsers[format](path)
#         if parser.check(rows=50):
#             return format
#     raise UnknownFormatError("Failed to infer the format, please specify it manually!")