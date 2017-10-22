"""Module for peaks."""

from bisect import bisect_left, bisect_right
from math import log, exp
from scipy.misc import comb
import random
import numpy as np
from statsmodels import api as sm


class Peak(object):
    """Class for peaks."""
    def __init__(self, c, s, e, smt=None):
        self.chrm = c
        self.start = s
        self.end = e
        self.summit = (s + e) / 2 + 1 if smt is None else smt + s
        self.read_count1 = 0.0
        self.read_density1 = 0.0
        self.normed_read_density1 = 0.0
        self.normed_read_density2 = 0.0
        self.read_count2 = 0.0
        self.read_density2 = 0.0
        self.mvalue = 0.0
        self.avalue = 0.0
        self.normed_mvalue = 0.0
        self.normed_avalue = 0.0
        self.pvalue = 0.0

    def set_summit(self, smt):
        self.summit = smt

    def __cal_read_count(self, reads_pos, ext):
        """Calculate read count of peak."""
        if self.chrm not in reads_pos.keys():
            return 0
        re_start, re_end = self.summit - ext - 1, self.summit + ext
        si = bisect_left(reads_pos[self.chrm], re_start)
        ei = bisect_right(reads_pos[self.chrm], re_end)
        try:
            if re_end == reads_pos[self.chrm][ei]:
                return ei - si + 1
            else:
                return ei - si
        except IndexError:
            return ei - si

    def __cal_read_density(self, reads_pos, ext):
        read_count = self.__cal_read_count(reads_pos, ext) + 1
        read_density = read_count * 1000.0 / (2.0 * ext)
        return read_count, read_density

    def cal_read_density(self, reads_pos1, reads_pos2, ext):
        self.read_count1, self.read_density1 = self.__cal_read_density(reads_pos1, ext)
        self.read_count2, self.read_density2 = self.__cal_read_density(reads_pos2, ext)
        self.mvalue = log(self.read_density1, 2) - log(self.read_density2, 2)
        self.avalue = (log(self.read_density1, 2) + log(self.read_density2, 2)) / 2

    def normalize_mavalue(self, ma_fit):
        """
        ma_fit: R2 = ma_fit[0] * R1 + ma_fit[1]
        """
        # key method for normalizing read density

        self.normed_mvalue = self.mvalue - (ma_fit[0] + ma_fit[1] * self.avalue)
        self.normed_avalue = self.avalue
        self.normed_read_density1 = 2 ** (self.normed_avalue + self.normed_mvalue / 2)
        self.normed_read_density2 = 2 ** (self.normed_avalue - self.normed_mvalue / 2)
        self.pvalue = _digit_exprs_p_norm(self.normed_read_density1, self.normed_read_density2)

    def isoverlap(self, other_pk):
        if self.start <= other_pk.start < self.end or self.start < other_pk.end <= self.end:
            return True
        else:
            return False


def _digit_exprs_p_norm(x, y):
    """Calculate p-value with p-value.
    """
    xx = round(x)
    if xx == 0:
        xx = 1
    yy = int(round(y))
    if yy == 0:
        yy = 1
    if xx + yy < 20.0:  # if x + y small
        p1 = round(comb(xx + yy, xx)) * 2 ** - (xx + yy + 1.0)
        p2 = round(comb(xx + yy, yy)) * 2 ** - (xx + yy + 1.0)
        return max(p1, p2)
    else:  # if x + y large, use the approximate equations
        log_p = (xx + yy) * log(xx + yy) - xx * log(xx) - yy * log(yy) - (xx + yy + 1.0) * log(2.0)
        if log_p < -500:
            log_p = -500
        p = exp(log_p)
        return p


def get_peaks_size(peaks):
    num = 0
    for chrom in peaks.keys():
        num += len(peaks[chrom])
    return num


def cal_peaks_read_density(peaks, reads_pos1, reads_pos2, ext):
    for chrom in peaks.keys():
        [peak.cal_read_density(reads_pos1, reads_pos2, ext) for peak in peaks[chrom]]


def normalize_peaks(peaks, ma_fit):
    for chrom in peaks.keys():
        [peak.normalize_mavalue(ma_fit) for peak in peaks[chrom]]


def get_common_peaks(peaks1, peaks2):
    peaks1_unique, peaks1_common, peaks2_unique, peaks2_common = {}, {}, {}, {}
    common_chrom = set(peaks1.keys()).intersection(peaks2.keys())
    peaks1_unique_chrom = set(peaks1.keys()).difference(common_chrom)
    peaks2_unique_chrom = set(peaks2.keys()).difference(common_chrom)
    for chrom in peaks1_unique_chrom:
        peaks1_unique[chrom] = peaks1[chrom]
    for chrom in peaks2_unique_chrom:
        peaks2_unique[chrom] = peaks2[chrom]
    for chrom in common_chrom:
        peaks1_unique[chrom], peaks1_common[chrom], peaks2_unique[chrom], peaks2_common[chrom] = \
            __get_common_peaks(peaks1[chrom], peaks2[chrom])
    return peaks1_unique, peaks1_common, peaks2_unique, peaks2_common


def __get_common_peaks(peaks1, peaks2):
    flag1, flag2 = np.zeros(len(peaks1)), np.zeros(len(peaks2))
    peaks2_start = np.array([peak.start for peak in peaks2])
    peaks2_end = np.array([peak.end for peak in peaks2])
    for idx, peak in enumerate(peaks1):
        claus = 1.0 * (peak.end - peaks2_start) * (peaks2_end - peak.start)
        overlap_locs = np.where(claus > 0)[0]
        if overlap_locs.size > 0:
            flag1[idx] = 1
            flag2[overlap_locs] = 1

    peaks1_unique, peaks1_common = [], []
    for idx, tmp_flag in enumerate(flag1):
        peaks1_unique.append(peaks1[idx]) if tmp_flag == 0 else peaks1_common.append(peaks1[idx])

    peaks2_unique, peaks2_common = [], []
    for idx, tmp_flag in enumerate(flag2):
        peaks2_unique.append(peaks2[idx]) if tmp_flag == 0 else peaks2_common.append(peaks2[idx])

    return peaks1_unique, peaks1_common, peaks2_unique, peaks2_common


def randomize_peaks(peaks):
    peaks_random = {}
    for chrom in peaks.keys():
        peaks_random[chrom] = []
        peaks_chrom = peaks[chrom]
        starts, ends = [peak.start for peak in peaks_chrom], [peak.end for peak in peaks_chrom]
        lengths = [e - s for e, s in zip(ends, starts)]
        min_start, max_end = min(starts), max(ends)
        for length in lengths:
            randomized_start = random.randint(min_start, max_end)
            peaks_random[chrom].append(Peak(chrom, randomized_start, randomized_start + length))
    return peaks_random


def merge_common_peaks(peaks1_common, peaks2_common):
    merged_pks = {}
    summit_dist = {}
    for chrom in peaks1_common.keys():
        mixed_peaks = peaks1_common[chrom] + peaks2_common[chrom]
        merged_pks[chrom], summit_dist[chrom] = __merge_sorted_peaks_list(_sort_peaks_list(mixed_peaks))
    return merged_pks, summit_dist


def _sort_peaks_list(peaks_list, by='start'):
    if by == 'start':
        return [peaks_list[loc] for loc in np.argsort([peak.start for peak in peaks_list])]
    elif by == 'summit':
        return [peaks_list[loc] for loc in np.argsort([peak.summit for peak in peaks_list])]


def _add_peaks(peaks1, peaks2):
    peaks = {}
    keys = set(peaks1.keys() + peaks2.keys())
    for chrom in keys:
        value = []
        if chrom in peaks1.keys():
            value += peaks1[chrom]
        if chrom in peaks2.keys():
            value += peaks2[chrom]
        peaks[chrom] = value
    return peaks


def __merge_sorted_peaks_list(sorted_peaks_list):
    # print 'sorted peaks length: %d' % len(sorted_pks_list)
    merged_peaks, smt_dists = [], []

    def get_a_merged_peak(head):
        # print head_loc
        merged_peaks_num = 0
        merged_peak = sorted_peaks_list[head]
        summits = [merged_peak.summit]
        for peak in sorted_peaks_list[head + 1:]:
            if merged_peak.isoverlap(peak):
                merged_peak = Peak(peak.chrm, min(merged_peak.start, peak.start), max(merged_peak.end, peak.end))
                summits.append(peak.summit)
                merged_peaks_num += 1
            else:
                sorted_summits = sorted(summits)
                smt_a, smt_b = get_summit(sorted_summits)
                smt_dist = smt_b - smt_a
                merged_peak.set_summit((smt_a + smt_b) / 2 + 1)
                new_head_loc = head + merged_peaks_num + 1
                return new_head_loc, merged_peak, smt_dist
        sorted_summits = sorted(summits)
        smt_a, smt_b = get_summit(sorted_summits)
        smt_dist = smt_b - smt_a
        merged_peak.set_summit((smt_a + smt_b) / 2 + 1)
        new_head_loc = head + merged_peaks_num + 1
        return new_head_loc, merged_peak, smt_dist

    h_loc = 0
    while h_loc < len(sorted_peaks_list):
        h_loc, m_pk, smt_d = get_a_merged_peak(h_loc)
        merged_peaks.append(m_pk)
        smt_dists.append(smt_d)
    return merged_peaks, smt_dists


def get_summit(sorted_summits):
    smt_starts, smt_ends = sorted_summits[:-1], sorted_summits[1:]
    smt_a, smt_b = smt_starts[0], smt_ends[0]
    for s, e in zip(smt_starts, smt_ends):
        if s - e < smt_a - smt_b:
            smt_a, smt_b = s, e
    return smt_a, smt_b


def use_merged_peaks_fit_model(merged_peaks, summit_dist, min_summit_dist):
    selected_pks = {}
    for chrom in merged_peaks.keys():
        selected_pks[chrom] = []
        for peaks, smt_d in zip(merged_peaks[chrom], summit_dist[chrom]):
            if smt_d <= min_summit_dist:
                selected_pks[chrom].append(peaks)
    mvalues, avalues = get_peaks_mavalues(selected_pks)
    fit_x = np.array(avalues)
    fit_y = np.array(mvalues)
    idx_sel = np.where((fit_y >= -10) & (fit_y <= 10))[0]
    # fit the model
    x = sm.add_constant(fit_x[idx_sel])
    y = fit_y[idx_sel]
    ma_fit = sm.RLM(y, x).fit().params
    return ma_fit


def get_peaks_mavalues(peaks):
    mvalues, avalues = [], []
    for chrom in peaks.keys():
        for peak in peaks[chrom]:
            mvalues.append(peak.mvalue), avalues.append(peak.avalue)
    return mvalues, avalues


def get_peaks_normed_mavalues(peaks):
    normed_mvalues, normed_avalues = [], []
    for chrom in peaks.keys():
        for peak in peaks[chrom]:
            normed_mvalues.append(peak.normed_mvalue), normed_avalues.append(peak.normed_avalue)
    return normed_mvalues, normed_avalues


def get_peaks_pvalues(peaks):
    pvalues = []
    for chrom in peaks.keys():
        for peak in peaks[chrom]:
            pvalues.append(peak.pvalue)
    return pvalues
