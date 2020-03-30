"""
manorm.region.utils
-------------------

Support functions for operations on MAnorm peaks.
"""

import logging
import random

import numpy as np

from manorm.region import GenomicRegion, ManormPeak, GenomicRegions

logger = logging.getLogger(__name__)


def overlap_on_same_chrom(regions1, regions2):
    """Given two sets of genomic regions(peaks) located on the same chromosome,
    returns the region overlap indicators of them.
    """
    overlap_flag1 = np.zeros(len(regions1), dtype=bool)
    overlap_flag2 = np.zeros(len(regions2), dtype=bool)
    for i, region_i in enumerate(regions1):
        for j, region_j in enumerate(regions2):
            if (region_i.end - region_j.start) * (
                    region_j.end - region_i.start) > 0:
                overlap_flag1[i] = True
                overlap_flag2[j] = True
    return overlap_flag1, overlap_flag2


def classify_peaks_by_overlap(peaks1, peaks2):
    """Classify two sets of peaks based on overlap and set the `iscommon` flag
    for every individual peak.
    """
    logger.debug("Classifying unique/common peaks by overlap")
    for chrom in set(peaks1.chroms) | set(peaks2.chroms):
        logger.debug(f"Classifying peaks on {chrom}")
        peaks1_chrom = peaks1.fetch(chrom)
        peaks2_chrom = peaks2.fetch(chrom)
        overlap_flag1, overlap_flag2 = overlap_on_same_chrom(
            peaks1_chrom, peaks2_chrom)
        for peak, flag in zip(peaks1_chrom, overlap_flag1):
            peak.iscommon = flag
        for peak, flag in zip(peaks2_chrom, overlap_flag2):
            peak.iscommon = flag
    return peaks1, peaks2


def merge_common_peaks(peaks1, peaks2):
    """Merge common (overlapping) peaks of the specified peak sets and
    returns the merged peaks.
    """

    def _merge_peak(chrom, start, end, summits):
        # take the middle of two nearest neighbour peaks as new summit
        summits.sort()
        min_dis = None
        for i, head in enumerate(summits[:-1]):
            tail = summits[i + 1]
            if min_dis is None or min_dis > (tail - head):
                min_dis = tail - head
                summit = (head + tail) // 2
        peak = ManormPeak(chrom, start, end, summit)
        peak.iscommon = True
        peak.summit_dis = min_dis
        return peak

    logger.debug("Merging common peaks")
    peaks_merged = GenomicRegions(name='merged_common_peaks')
    for chrom in set(peaks1.chroms) & set(peaks2.chroms):
        logger.debug(f"Merging peaks on {chrom}")
        peaks_mixed = []
        for peak in peaks1.fetch(chrom):
            if peak.iscommon:
                peaks_mixed.append(peak)
        for peak in peaks2.fetch(chrom):
            if peak.iscommon:
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
                    peaks_merged.add(_merge_peak(chrom, start, end, summits))
                    start = peak.start
                    end = peak.end
                    summits = [peak.summit]
            else:  # add the final one
                peaks_merged.add(_merge_peak(chrom, start, end, summits))
                break
    return peaks_merged


def generate_random_regions(ref_regions):
    """Generate random control regions from the given reference regions.
    The length and chromosome distribution for each region are controlled to
    match the reference regions.
    """
    regions_random = GenomicRegions(name='random')
    for chrom in ref_regions.chroms:
        templates = ref_regions.fetch(chrom)
        if len(templates) == 0:
            continue
        starts = [region.start for region in templates]
        lengths = [region.end - region.start for region in templates]
        start_min = min(starts)
        start_max = max(starts)
        for length in lengths:
            tmp_start = random.randint(start_min, start_max)
            regions_random.add(
                GenomicRegion(chrom, tmp_start, tmp_start + length))
    return regions_random


def random_peak_overlap(peaks1, peaks2, n_random):
    """Calculate the number of overlapping peaks between peaks1 and
    random control peaks generated based on peaks2.
    """
    n_overlap_rand = []
    for _ in range(n_random):
        peak_rand = generate_random_regions(peaks2)
        n_overlap = 0
        for chrom in set(peaks1.chroms) & set(peak_rand.chroms):
            flag_overlap, _ = overlap_on_same_chrom(peaks1.fetch(chrom),
                                                    peak_rand.fetch(chrom))
            n_overlap += flag_overlap.sum()
        n_overlap_rand.append(n_overlap)
    n_overlap_rand = np.array(n_overlap_rand)
    n_overlap_mean = n_overlap_rand.mean()
    n_overlap_std = n_overlap_rand.std()
    return n_overlap_mean, n_overlap_std


def count_common_peaks(peaks):
    """Returns the number of common peaks."""
    n_common = 0
    for chrom in peaks.chroms:
        for peak in peaks.fetch(chrom):
            if peak.iscommon:
                n_common += 1
    return n_common


def count_unique_peaks(peaks):
    """Returns the number of unique peaks."""
    n_unique = 0
    for chrom in peaks.chroms:
        for peak in peaks.fetch(chrom):
            if not peak.iscommon:
                n_unique += 1
    return n_unique
