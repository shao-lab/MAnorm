# -*- coding: utf-8 -*-

"""
manorm.read
~~~~~~~~~~~

This module contains Reads objects for read-related operations.
"""

import os
import logging
from bisect import bisect_left

from manorm.exceptions import UnmatchedBedFormatError, UnsupportedFormatError
from manorm.read.parsers import BAMParser, BEDPEParser, BEDParser, SAMParser

READ_FORMATS = ['bed', 'bedpe', 'sam', 'bam']
READ_PARSERS = {'bed': BEDParser, 'bedpe': BEDPEParser, 'sam': SAMParser, 'bam': BAMParser}

logger = logging.getLogger(__name__)


class Reads(object):
    """Class for reads generated from next-generation sequencing."""

    def __init__(self, name=None):
        """Initialize the read set.

        :param name: The name of the read set.
        """
        self.name = name
        self.data = {}

    @property
    def chroms(self):
        """Return the chromosome names of reads."""
        return self.data.keys()

    @property
    def size(self):
        """Return the number of reads."""
        return sum(len(self.data[chrom]) for chrom in self.chroms)

    def add(self, chrom, pos):
        """Add a read.

        :param chrom: Chromosome name of the read.
        :param pos: Representative position of the read.
        """
        self.data.setdefault(chrom, [])
        self.data[chrom].append(pos)

    def sort(self):
        """Sort reads."""
        for chrom in self.chroms:
            self.data[chrom].sort()

    def count(self, chrom, start, end):
        """Count reads in given interval by binary search."""
        try:
            head = bisect_left(self.data[chrom], start)
            tail = bisect_left(self.data[chrom], end)
            return tail - head
        except KeyError:
            return 0


def load_reads(path, format='bed', paired=False, shift=100, name=None):
    """Read reads from file.

    :param path: The file path to read reads from.
    :param format: Format of reads file.
    :param paired: Paired-end mode.
    :param shift: Shift size of single-end reads.
    :param name: Name of reads.
    """
    logger.debug("Loading reads from {}".format(path))
    if paired and format == 'bed':
        raise UnmatchedBedFormatError
    if not paired and format == 'bedpe':
        raise UnmatchedBedFormatError
    if name is None:
        name = os.path.splitext(os.path.basename(path))[0]
    reads = Reads(name=name)
    try:
        read_parser = READ_PARSERS[format](path)
    except KeyError:
        raise UnsupportedFormatError(format=format)
    for chrom, pos in read_parser.parse(paired=paired, shift=shift):
        reads.add(chrom, pos)
    read_parser.close()
    reads.sort()
    logger.debug("Loaded {:,} reads".format(reads.size))
    return reads
