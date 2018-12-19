# -*- coding: utf-8 -*-

"""
manorm.peak.utils
~~~~~~~~~~~~~~~~~

Support functions for peaks.
"""

from __future__ import absolute_import, division

import logging
import random

import numpy as np

from manorm.peak import Peak, Peaks

logger = logging.getLogger(__name__)


def overlap_on_single_chr(peaks1, peaks2):
    """Given two sets of peaks on the same chromosome, return the overlap flag."""
    flag1 = np.zeros(len(peaks1), dtype=np.int)
    flag2 = np.zeros(len(peaks2), dtype=np.int)
    starts = np.array([peak.start for peak in peaks2])
    ends = np.array([peak.end for peak in peaks2])
    for idx, peak in enumerate(peaks1):
        product = (peak.end - starts) * (ends - peak.start)
        overlap_idx = np.where(product > 0)[0]
        if overlap_idx.size > 0:
            flag1[idx] = 1
            flag2[overlap_idx] = 1
    return flag1, flag2


def classify_peaks_by_overlap(peaks1, peaks2):
    """Given two sets of peaks, classify peaks according to overlap."""
    logger.debug("Classifying unique/common peaks by overlap")
    for chrom in set(peaks1.chroms) | set(peaks2.chroms):
        logger.debug("Classifying on {}".format(chrom))
        flag1, flag2 = overlap_on_single_chr(peaks1.fetch(chrom), peaks2.fetch(chrom))
        for idx, temp_flag in enumerate(flag1):
            if temp_flag == 0:
                peaks1.data[chrom][idx].type = 'unique'
            else:
                peaks1.data[chrom][idx].type = 'common'
        for idx, temp_flag in enumerate(flag2):
            if temp_flag == 0:
                peaks2.data[chrom][idx].type = 'unique'
            else:
                peaks2.data[chrom][idx].type = 'common'
    return peaks1, peaks2


def merge_common_peaks(peaks1, peaks2):
    """Merge common peaks of two peak sets and return merged peaks."""

    def _add_merged_peak(chrom, start, end, summits):
        min_dis = None
        summits.sort()
        # the summit of merged peaks is the middle point of two nearest neighbour peaks
        for i, head in enumerate(summits[:-1]):
            tail = summits[i + 1]
            if min_dis is None or tail - head < min_dis:
                min_dis = tail - head
                summit = (head + tail) // 2
        peak = Peak(chrom, start, end, summit)
        peak.type = "merged_common_peaks"
        peak.summit_dis = min_dis
        peaks_merged.add(peak)

    logger.debug("Merging common peaks")
    peaks_merged = Peaks(name='merged_common_peaks')

    for chrom in set(peaks1.chroms) & set(peaks2.chroms):
        logger.debug("Merging on {}".format(chrom))
        # mix common peaks together
        peaks_mixed = []
        for peak in peaks1.fetch(chrom):
            if peak.type == 'common':
                peaks_mixed.append(peak)
        for peak in peaks2.fetch(chrom):
            if peak.type == 'common':
                peaks_mixed.append(peak)
        if len(peaks_mixed) == 0:
            continue

        peaks_mixed.sort(key=lambda x: x.start)

        idx = 0
        peak = peaks_mixed[idx]
        start = peak.start  # the left-most pos of continuous overlapping peaks
        end = peak.end  # the right-most pos of continuous overlapping peaks
        summits = [peak.summit]
        while True:
            idx += 1
            if idx < len(peaks_mixed):
                peak = peaks_mixed[idx]
                if peak.start < end:  # current peak overlaps with the last one
                    end = max(end, peak.end)
                    summits.append(peak.summit)
                else:  # not overlap
                    _add_merged_peak(chrom, start, end, summits)
                    start = peak.start
                    end = peak.end
                    summits = [peak.summit]
            else:  # add the final merged peak
                _add_merged_peak(chrom, start, end, summits)
                break
    return peaks_merged


def generate_random_peaks(ref):
    """Generate a set of random control peaks from given reference peaks.
    The random peaks have matched peak lengths and peak counts for each chromosome with reference peaks.
    """
    peaks_rand = Peaks(name='random')
    for chrom in ref.chroms:
        temp_peaks = ref.fetch(chrom)
        if len(temp_peaks) == 0:
            continue
        starts = [peak.start for peak in temp_peaks]
        lengths = [peak.end - peak.start for peak in temp_peaks]
        start_min = min(starts)
        start_max = max(starts)
        for length in lengths:
            start_rand = random.randint(start_min, start_max)
            peaks_rand.add(Peak(chrom, start_rand, start_rand + length))
    return peaks_rand
