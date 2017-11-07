"""Module for reads and related operation."""
import os
from collections import defaultdict
from manorm.logger import logger


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
            self.load_reads(path, shift)
            self.sort()

    def load_reads(self, path, shift):
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
                self.pos[chrom].append(pos)
        logger.info("{} reads are loaded from {}".format(self.size, os.path.basename(path)))

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
