"""Module for peaks and related operation."""
import os
import random
from math import log, exp
import numpy as np
from bisect import bisect_left
from collections import defaultdict
from manorm.logger import logger


class Peaks(object):
    """Class for a set of peaks."""

    def __init__(self, path=None, name=None):
        self.path = path
        self.name = name
        if self.path and not self.name:
            self.name = os.path.splitext(os.path.basename(path))[0]

        self.peaks = defaultdict(list)
        self._size = None
        if self.path:
            self.load_peaks(path)

    def load_peaks(self, path):
        """Load peaks."""

        def _bed_parser(path):
            """Parser for peaks."""
            with open(path, "r") as fin:
                for line in fin:
                    if line.startswith("#"):
                        continue
                    fields = line.split("\t")
                    chrom = fields[0].strip()
                    try:
                        peak = Peak(chrom, int(fields[1]), int(fields[2]), int(fields[3]))
                    except:
                        peak = Peak(chrom, int(fields[1]), int(fields[2]))
                    self.peaks[chrom].append(peak)

        def _macs_xls_parser(path):
            """Parser for macs xls peaks."""
            with open(path, "r") as fin:
                for line in fin:
                    if line.startswith("#"):
                        continue
                    fields = line.split("\t")
                    chrom = fields[0]
                    try:
                        peak = Peak(chrom, int(fields[1]) - 1, int(fields[2]), int(fields[4]))
                    except:
                        continue
                    self.peaks[chrom].append(peak)

        try:
            _bed_parser(path)
        except:
            self.peaks = defaultdict(list)
            _macs_xls_parser(path)
        logger.info("{} peaks are loaded from {}".format(self.size, os.path.basename(path)))

    @property
    def size(self):
        if not self._size:
            num = 0
            for chrom in self.peaks:
                num += len(self.peaks[chrom])
            self._size = num
        return self._size

    def count_peak_type(self, key):
        num = 0
        for chrom in self.peaks:
            for peak in self.peaks[chrom]:
                if peak.type and peak.type.endswith(key):
                    num += 1
        return num


class Peak(object):
    """Class for peak."""

    def __init__(self, chrom, start, end, summit=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        if summit:
            if start + summit < end:  # relative position
                self.summit = start + summit
            else:  # absolute position
                self.summit = summit
        else:
            self.summit = (start + end) / 2
        self.type = None
        self.read_count1 = None
        self.read_count2 = None
        self.read_density1 = None
        self.read_density2 = None
        self.normed_read_density1 = None
        self.normed_read_density2 = None
        self.m_value = None
        self.a_value = None
        self.normed_m_value = None
        self.normed_a_value = None
        self.p_value = None

    def _count_reads(self, reads, extend):
        """Count reads by binary search."""
        try:
            head = bisect_left(reads.pos[self.chrom], self.summit - extend)
            tail = bisect_left(reads.pos[self.chrom], self.summit + extend)
            return tail - head
        except KeyError:
            return 0

    def cal_read_count(self, reads1, reads2, extend):
        """Calculate the read count and read density."""
        self.read_count1 = self._count_reads(reads1, extend) + 1
        self.read_count2 = self._count_reads(reads2, extend) + 1
        self.read_density1 = self.read_count1 * 1000.0 / (extend * 2)
        self.read_density2 = self.read_count2 * 1000.0 / (extend * 2)

    def cal_ma_value(self):
        """Calculate the M value and A value."""
        self.m_value = log(self.read_density1, 2) - log(self.read_density2, 2)
        self.a_value = (log(self.read_density1, 2) + log(self.read_density2, 2)) / 2

    def normalize(self, ma_model):
        """Normalize M value and A value by robust linear model.
        ma_model: y = ma_model[0] * x + ma_model[1]
        """

        def _cal_p_value(x, y):
            """Calculate P value with given read densities."""

            def _combination(n, k):
                import operator as op
                if n < 0 or k < 0 or k > n:
                    value = 0
                else:
                    k = min(k, n - k)
                    if k == 0:
                        value = 1
                    else:
                        numerator = reduce(op.mul, range(n, n - k, -1))
                        denominator = reduce(op.mul, range(1, k + 1))
                        value = numerator // denominator
                return value

            def _log_factorial(n):
                num = 0
                for i in range(1, n + 1):
                    num += log(i)
                return num

            x = int(round(x))
            if x == 0:
                x = 1
            y = int(round(y))
            if y == 0:
                y = 1
            if x + y < 20:
                p_value = _combination(x + y, x) * 2 ** - (x + y + 1)
            else:  # if x + y is large, use the log-transform to calculate p-value
                log_p = _log_factorial(x + y) - _log_factorial(x) - _log_factorial(y) - (x + y + 1) * log(2)
                # log_p = (x + y) * log(x + y) - x * log(x) - y * log(y) - (x + y + 1) * log(2)
                if log_p < -500:
                    log_p = -500
                p_value = exp(log_p)
            return p_value

        self.normed_m_value = round(self.m_value - (ma_model[0] + ma_model[1] * self.a_value), 5)
        self.normed_a_value = round(self.a_value, 5)
        self.normed_read_density1 = round(2 ** (self.normed_a_value + self.normed_m_value / 2), 5)
        self.normed_read_density2 = round(2 ** (self.normed_a_value - self.normed_m_value / 2), 5)
        self.p_value = _cal_p_value(self.normed_read_density1, self.normed_read_density2)


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
        merged_peaks.peaks[chrom] = temp_merged_peaks
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
            random_peaks.peaks[chrom].append(Peak(chrom, random_start, random_start + length))
    return random_peaks
