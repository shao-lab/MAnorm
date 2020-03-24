"""
manorm.read
-----------

Module for sequencing reads.
"""

import logging
import os
from bisect import bisect_left

from manorm.exceptions import FormatModeConflictError
from manorm.read.parsers import get_read_parser

READ_FORMATS = ['bed', 'bedpe', 'sam', 'bam']

logger = logging.getLogger(__name__)


class Reads:
    """Class for reads generated from next-generation sequencing.

    Parameters
    ----------
    name : str, optional
        The sample name of the sequencing reads.

    Attributes
    ----------
    name : str or None
        The sample name of the sequencing reads.
    """

    def __init__(self, name=None):
        self.name = name
        self._data = {}

    @property
    def chroms(self):
        """Returns sorted chromosome names of the sequencing reads.

        Returns
        -------
        list of str
            Chromosome names (sorted) of the sequencing reads.
        """
        return sorted(self._data.keys())

    @property
    def size(self):
        """Returns the number of sequencing reads."""
        return sum(len(value) for value in self._data.values())

    def add(self, chrom, pos):
        """Add a read position.

        Parameters
        ----------
        chrom : str
            The chromosome name of the read.
        pos : int
            The representative genomic position of the read.
        """
        self._data.setdefault(chrom, [])
        self._data[chrom].append(pos)

    def sort(self):
        """Sort reads."""
        for chrom in self.chroms:
            self._data[chrom].sort()

    def count(self, chrom, start, end):
        """Count reads located in the given interval by binary search.

        Parameters
        ----------
        chrom : str
            The chromosome name of the interval.
        start : int
            The start pos of the interval.
        end : int
            The end pos of the interval.
        """
        if start >= end:
            raise ValueError(
                f"expect start < end, got: start={start} end={end}")
        try:
            head = bisect_left(self._data[chrom], start)
            tail = bisect_left(self._data[chrom], end)
            return tail - head
        except KeyError:
            return 0


def load_reads(path, format='bed', paired=False, shift=100, name=None):
    """Read reads from file.

    Parameters
    ----------
    path : str
        Path to load the reads.
    format : str, optional
        File format, default='bed'.
    paired : bool, optional
        Whether the reads are paired-end or not, default=False.
    shift : int, optional
        Shift size for single-end reads, default=100.
    name : str, optional
        Sample name. If not specified, the basename of the file will be used.

    Returns
    -------
    reads : `Reads`
        Loaded sequencing reads.
    """
    logger.info(f"Loading reads from {path} [{format}]")
    if format == 'bed' and paired:
        raise FormatModeConflictError('bed', 'paired-end')
    if format == 'bedpe' and not paired:
        raise FormatModeConflictError('bedpe', 'single-end')
    if name is None:
        name = os.path.splitext(os.path.basename(path))[0]
    reads = Reads(name=name)
    parser = get_read_parser(format)(path)
    for chrom, pos in parser.parse(paired=paired, shift=shift):
        reads.add(chrom, pos)
    reads.sort()
    logger.info(f"Loaded {reads.size:,} reads")
    return reads
