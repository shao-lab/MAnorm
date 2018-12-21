# -*- coding: utf-8 -*-

"""
manorm.peak
~~~~~~~~~~~

This module contains classes and functions for peak-related operations.
"""

from __future__ import absolute_import, division

import logging
import os
from math import exp, log

from manorm.compat import filter, range
from manorm.exceptions import UnsupportedFormatError
from manorm.peak.parsers import BEDParser, BEDSummitParser, MACS2Parser, MACSParser, NarrowPeakParser

logger = logging.getLogger(__name__)

PEAK_FORMATS = ['bed', 'bed-summit', 'macs', 'macs2', 'narrowpeak', 'broadpeak']
PEAK_PARSERS = {'bed': BEDParser, 'bed-summit': BEDSummitParser, 'macs': MACSParser, 'macs2': MACS2Parser,
                'narrowpeak': NarrowPeakParser}


class Peak(object):
    """Class for a single peak."""

    def __init__(self, chrom, start, end, summit=None):
        """Initialize a peak.

        :param chrom: The chromosome name of the peak.
        :param start: The start coordinate of the peak (0-based).
        :param end: The end coordinate of the peak (0-based).
        :param summit: The summit coordinate of the peak (0-based).
        """
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        if summit is not None:
            self.summit = int(summit)
        else:
            self.summit = (self.start + self.end) // 2
        if not self.start <= self.summit <= self.end:
            raise ValueError("peak start must be <= summit and < end, got start={}, summit={}, end={}".format(
                self.start, self.summit, self.end))
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

    def count_reads(self, reads1, reads2, window_size=2000):
        """Count reads and calculate the read density."""
        if window_size <= 0:
            raise ValueError("window size must be > 0")
        extend = window_size // 2
        self.read_count1 = reads1.count(self.chrom, self.summit - extend, self.summit + extend) + 1  # add a pseudo 1
        self.read_count2 = reads2.count(self.chrom, self.summit - extend, self.summit + extend) + 1
        self.read_density1 = self.read_count1 * 1000 / (extend * 2)
        self.read_density2 = self.read_count2 * 1000 / (extend * 2)

    def cal_ma_value(self):
        """Calculate the M value and A value based on read densities."""
        self.m_value = log(self.read_density1, 2) - log(self.read_density2, 2)
        self.a_value = (log(self.read_density1, 2) + log(self.read_density2, 2)) / 2

    def normalize(self, ma_params):
        """Normalize M value and A value by robust linear model.
        ma_model: y = ma_params[1] * x + ma_params[0]
        """

        def _cal_p_value(x, y):
            """Calculate P value with given read densities."""

            def _log_factorial(n):
                num = 0
                for i in range(1, n + 1):
                    num += log(i)
                return num

            if x < 0 or y < 0:
                raise ValueError("x and y must be >= 0")
            x = int(round(x))
            if x == 0:
                x = 1
            y = int(round(y))
            if y == 0:
                y = 1
            # use the log-transform to calculate p-value
            log_p = _log_factorial(x + y) - _log_factorial(x) - _log_factorial(y) - (x + y + 1) * log(2)
            if log_p < -500:
                log_p = -500
            p_value = exp(log_p)
            return p_value

        self.m_value_normed = round(self.m_value - (ma_params[0] + ma_params[1] * self.a_value), 5)
        self.a_value_normed = round(self.a_value, 5)
        self.read_density1_normed = round(2 ** (self.a_value_normed + self.m_value_normed / 2), 5)
        self.read_density2_normed = round(2 ** (self.a_value_normed - self.m_value_normed / 2), 5)
        self.p_value = _cal_p_value(self.read_density1_normed, self.read_density2_normed)
        self.normed = True

    def __repr__(self):
        return "Peak({}:{}-{})".format(self.chrom, self.start, self.end)


class Peaks(object):
    """Class for a collection of peaks.
    Peaks are stored by a dict under ``self.data`` with chromosome names as keys and lists of :class:`Peak` as values.
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
        return list(self.data.keys())

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

        :param by: Attribute name to sort by. Defaults to ``'start'``.
        :param ascending: Sort ascending or descending. Defaults to ``True``.
        """
        for chrom in self.chroms:
            self.data[chrom].sort(key=lambda x: getattr(x, by), reverse=not ascending)

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
        raise UnsupportedFormatError(format=format)
    for chrom, start, end, summit in peak_parser.parse():
        peaks.add(Peak(chrom=chrom, start=start, end=end, summit=summit))
    peak_parser.close()
    peaks.sort()
    logger.debug("Loaded {} peaks".format(peaks.size))
    return peaks
