# -*- coding: utf-8 -*-

"""
manorm.peak
~~~~~~~~~~~

This module contains Peak and Peaks objects for peak-related operations.
"""

from __future__ import absolute_import, division

import logging
import os

from manorm.compat import filter
from manorm.exceptions import UnknownFormatError
from manorm.peak.parsers import BEDParser, BEDSummitParser, MACS2Parser, MACSParser, NarrowPeakParser

logger = logging.getLogger(__name__)

PEAK_FORMATS = ['bed', 'bed-summit', 'macs', 'macs2', 'narrowpeak', 'broadpeak']
PEAK_PARSERS = {'bed': BEDParser, 'bed-summit': BEDSummitParser, 'macs': MACSParser, 'macs2': MACS2Parser,
                'narrowpeak': NarrowPeakParser}


class Peak(object):
    """Class for a single peak."""

    def __init__(self, chrom, start, end, summit=None):
        """Initialize a peak.

        :param chrom: Chromosome name.
        :param start: Start coordinate (0-based).
        :param end: End coordinate (0-based).
        :param summit: Summit coordinate (0-based).
        """
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        if summit is not None:
            self.summit = int(summit)
        else:
            self.summit = (self.start + self.end) // 2
        self.type = None
        self.summit_dis = None
        self.read_count1 = None
        self.read_count2 = None
        self.read_density1 = None
        self.read_density2 = None
        self.m_value = None
        self.a_value = None
        self.normed = False
        self.read_density1_normed = None
        self.read_density2_normed = None
        self.m_value_normed = None
        self.a_value_normed = None
        self.p_value = None

    def __repr__(self):
        return "Peak({}:{}-{})".format(self.chrom, self.start, self.end)


class Peaks(object):
    """Class for a collection of peaks.
    Peaks are stored by a dict under `self.data` with chromosome names as keys and lists of :class:`Peak` as values.
    """

    def __init__(self, name=None):
        """Initialize the peak set.

        :param name: The name of the peak set.
        """
        self.name = name
        self.data = {}

    @property
    def chroms(self):
        """Return the chromosome names of peaks."""
        return self.data.keys()

    @property
    def size(self):
        """Return the number of peaks."""
        return sum(len(self.data[chrom]) for chrom in self.chroms)

    def add(self, peak):
        """Add a peak.

        :param peak: An instance of :class：`Peak` to be added.
        """
        if not isinstance(peak, Peak):
            raise ValueError("requires a 'Peak' object to be added into peaks")
        else:
            self.data.setdefault(peak.chrom, [])
            self.data[peak.chrom].append(peak)

    def sort(self, by='start', ascending=True):
        """Sort peaks.

        :param by: Attribute name to sort by. Defaults to `'start'`.
        :param ascending: Sort ascending or descending. Defaults to `True`.
        """
        for chrom in self.chroms:
            self.data[chrom].sort(key=lambda x: getattr(x, by), reverse=ascending)

    def fetch(self, chrom):
        """Fetch peaks from specified chromosome.

        :param chrom: Chromosome name to fetch peaks from.
        """
        if chrom in self.data:
            return self.data[chrom]
        else:
            return []

    @property
    def n_common(self):
        """Return the number of common peaks."""
        return sum(sum(1 for _ in filter(lambda x: x.type == 'common', self.fetch(chrom))) for chrom in self.chroms)

    @property
    def n_unique(self):
        """Return the number of unique peaks."""
        return sum(sum(1 for _ in filter(lambda x: x.type == 'unique', self.fetch(chrom))) for chrom in self.chroms)

    def __repr__(self):
        return "Peaks(name={})".format(self.name)


def load_peaks(path, format='bed', name=None):
    """Read peaks from file.

    :param path: The file path to read peaks from.
    :param format: Format of peaks file.
    :param name: Name of peaks.
    """
    logger.debug("Loading peaks from {}".format(path))
    if name is None:
        name = os.path.splitext(os.path.basename(path))[0]
    peaks = Peaks(name=name)
    try:
        peak_parser = PEAK_PARSERS[format](path)
    except KeyError:
        raise UnknownFormatError(format=format)
    for chrom, start, end, summit in peak_parser.parse():
        peaks.add(Peak(chrom=chrom, start=start, end=end, summit=summit))
    logger.debug("Loaded {} peaks".format(peaks.size))
    return peaks
