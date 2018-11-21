# -*- coding: utf-8 -*-

"""
manorm.peak
~~~~~~~~~~~

This module contains Peak and Peaks objects for peak-related operations.
"""

from __future__ import absolute_import, division

import logging
import random
import os
import numpy as np

logger = logging.getLogger(__name__)

PEAK_FORMATS = ['bed', 'macs', 'macs2', 'narrowpeak', 'broadpeak']


class Peak(object):
    """A single peak."""

    def __init__(self, chrom, start, end, summit=None):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        if summit is not None:
            self.summit = int(summit)
        else:
            self.summit = (self.start + self.end) // 2
        self.type = None
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
    """A collection of peaks.
    Peaks are stored into a dict with chromosome names as keys and lists of :class:`Peak` as values.
    """

    def __init__(self, name=None):
        self.name = name
        self.data = {}

    @property
    def size(self):
        """Return the number of peaks."""
        return sum(len(self.data[chrom]) for chrom in self.data)

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
        for chrom in self.data:
            self.data[chrom].sort(key=lambda x: getattr(x, by), reverse=ascending)

    def fetch(self, chrom):
        """Fetch peaks from specified chromosome.

        :param chrom: Chromosome name to fetch peaks from.
        """
        if chrom in self.data:
            return self.data[chrom]
        else:
            return []

    def __repr__(self):
        return "Peaks(name={})".format(self.name)


def load_peaks(path, format='auto', name=None):
    """Load peaks from file.

    :param path: The file path to read peaks from.
    :param format: Format of peaks file. Defaults to ``'auto'``.
    :param name: Name of peaks.
    """
    from manorm.peak.parsers import get_peak_parser
    if name is None:
        name = os.path.splitext(os.path.basename(path))[0]
    peaks = Peaks(name=name)
    peak_parser = get_peak_parser(path=path, format=format)
    for peak in peak_parser.parse():
        peaks.add(peak)
    return peaks


def overlap_on_single_chr(peaks1, peaks2):
    overlap_flag1 = np.zeros(len(peaks1))
    overlap_flag2 = np.zeros(len(peaks2))
    starts = np.array([peak.start for peak in peaks2])
    ends = np.array([peak.end for peak in peaks2])
    for idx, peak in enumerate(peaks1):
        product = (peak.end - starts) * (ends - peak.start)
        overlap_idx = np.where(product > 0)[0]
        if overlap_idx.size > 0:
            overlap_flag1[idx] = 1
            overlap_flag2[overlap_idx] = 1
    return overlap_flag1, overlap_flag2


def merge_peaks(peaks1, peaks2):
    def _merge_summits(summits):
        summits.sort()
        min_dis = None
        for idx, head in enumerate(summits[:-1]):
            tail = summits[idx + 1]
            if not min_dis or tail - head < min_dis:
                summit = (head + tail) / 2
                min_dis = tail - head
        return summit, min_dis

    merged_peaks = Peaks(name="merged_common_peaks")
    summit_dis = defaultdict(list)
    for chrom in set(peaks1) | set(peaks2):
        mixed_peaks = sorted(peaks1[chrom] + peaks2[chrom], key=lambda x: x.start)
        if len(mixed_peaks) == 0:
            logger.warning("Mixed peaks are empty on {}".format(chrom))
            continue
        temp_merged_peaks = []
        temp_summit_dis = []
        start = mixed_peaks[0].start
        end = mixed_peaks[0].end
        summits = [mixed_peaks[0].summit]
        idx = 1
        while idx < len(mixed_peaks):
            peak = mixed_peaks[idx]
            if peak.start < end:
                end = max(end, peak.end)
                summits.append(peak.summit)
            else:
                summit, dis = _merge_summits(summits)
                temp_merged_peaks.append(Peak(chrom, start, end, summit))
                temp_summit_dis.append(dis)
                start = peak.start
                end = peak.end
                summits = [peak.summit]
            idx += 1
        summit, dis = _merge_summits(summits)
        temp_merged_peaks.append(Peak(chrom, start, end, summit))
        temp_summit_dis.append(dis)
        for peak in temp_merged_peaks:
            peak.type = "merged_common_peaks"
        merged_peaks.data[chrom] = temp_merged_peaks
        summit_dis[chrom] = temp_summit_dis
    return merged_peaks, summit_dis


def generate_random_peaks(peaks):
    random_peaks = Peaks(name="random")
    for chrom in peaks:
        temp_peaks = peaks[chrom]
        if not temp_peaks:
            continue
        starts = [peak.start for peak in temp_peaks]
        lengths = [peak.end - peak.start for peak in temp_peaks]
        min_start = min(starts)
        max_start = max(starts)
        for length in lengths:
            random_start = random.randint(min_start, max_start)
            random_peaks.data[chrom].append(Peak(chrom, random_start, random_start + length))
    return random_peaks


def classify_peaks_by_overlap(peaks1, peaks2):
    for chrom in set(self.peaks1.peaks) | set(self.peaks2.peaks):
        overlap_flag1, overlap_flag2 = overlap_on_single_chr(self.peaks1.peaks[chrom], self.peaks2.peaks[chrom])
        for idx, temp_flag in enumerate(overlap_flag1):
            if temp_flag == 0:
                self.peaks1.peaks[chrom][idx].type = self.peaks1.name + "_unique"
            else:
                self.peaks1.peaks[chrom][idx].type = self.peaks1.name + "_common"
        for idx, temp_flag in enumerate(overlap_flag2):
            if temp_flag == 0:
                self.peaks2.peaks[chrom][idx].type = self.peaks2.name + "_unique"
            else:
                self.peaks2.peaks[chrom][idx].type = self.peaks2.name + "_common"
