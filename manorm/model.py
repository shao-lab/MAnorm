"""This module contains Peak and Peaks objects for peak-related operations."""

from __future__ import absolute_import
from manorm.peak import Peaks
from manorm.read import Reads

def cal_read_count(self, reads1, reads2, extend):
    """Calculate the read count and read density.

    :param reads1:
    :param reads2:
    :param extend:
    """
    self.read_count1 = reads1.count(self.chrom, self.summit -extend, self.summit +extend) + 1
    self.read_count2 = reads2.count(self.chrom, self.summit -extend, self.summit +extend) + 1
    self.read_density1 = self.read_count1 * 1000 / (extend * 2)
    self.read_density2 = self.read_count2 * 1000 / (extend * 2)

def cal_ma_value(self):
    """Calculate the M value and A value based on read densities."""
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

    self.m_value_normed = round(self.m_value - (ma_model[0] + ma_model[1] * self.a_value), 5)
    self.a_value_normed = round(self.a_value, 5)
    self.read_density1_normed = round(2 ** (self.a_value_normed + self.m_value_normed / 2), 5)
    self.read_density2_normed = round(2 ** (self.a_value_normed - self.m_value_normed / 2), 5)
    self.p_value = _cal_p_value(self.read_density1_normed, self.read_density2_normed)
    self.normed = True


class MAmodel(object):
    def __init__(self, peaks1, peaks2, reads1, reads2):
        self.peaks1 = peaks1
        self.peaks2 = peaks2
        self.reads1 = reads1
        self.reads2 = reads2


    def merge_common_peaks(self):
        peaks1_common = defaultdict(list)
        peaks2_common = defaultdict(list)
        for chrom in self.peaks1.peaks:
            for peak in self.peaks1.peaks[chrom]:
                if peak.type.endswith("common"):
                    peaks1_common[chrom].append(peak)
        for chrom in self.peaks2.peaks:
            for peak in self.peaks2.peaks[chrom]:
                if peak.type.endswith("common"):
                    peaks2_common[chrom].append(peak)
        self.common_peaks, self.common_dist = merge_peaks(peaks1_common, peaks2_common)

    def overlap_enrichment_test(self, n_random):
        real_overlap_num = self.peaks1.count_peak_type("common")
        random_overlap_num = []
        for _ in range(n_random):
            random_peaks = generate_random_peaks(self.peaks2.peaks)
            temp_overlap_num = 0
            for chrom in set(self.peaks1.peaks) | set(random_peaks.peaks):
                overlap_flag, _ = overlap_on_single_chr(self.peaks1.peaks[chrom], random_peaks.peaks[chrom])
                temp_overlap_num += overlap_flag.sum()
            if temp_overlap_num == 0:
                temp_overlap_num = 1.0
            random_overlap_num.append(temp_overlap_num)
        random_overlap_num = np.array(random_overlap_num)
        fold_change = real_overlap_num / random_overlap_num
        return random_overlap_num, fold_change

    def cal_read_density(self, extend=0):
        for chrom in self.peaks1.peaks:
            for peak in self.peaks1.peaks[chrom]:
                peak.cal_read_count(self.reads1, self.reads2, extend)
                peak.cal_ma_value()
        for chrom in self.peaks2.peaks:
            for peak in self.peaks2.peaks[chrom]:
                peak.cal_read_count(self.reads1, self.reads2, extend)
                peak.cal_ma_value()
        for chrom in self.common_peaks.peaks:
            for peak in self.common_peaks.peaks[chrom]:
                peak.cal_read_count(self.reads1, self.reads2, extend)
                peak.cal_ma_value()

    def fit_model(self, summit_dis_cutoff):
        m_values = []
        a_values = []
        for chrom in self.common_peaks.peaks:
            for idx, peak in enumerate(self.common_peaks.peaks[chrom]):
                if self.common_dist[chrom][idx] <= summit_dis_cutoff:
                    m_values.append(peak.m_value)
                    a_values.append(peak.a_value)
        m_values = np.array(m_values)
        a_values = np.array(a_values)
        mask = np.where((m_values >= -10) & (m_values <= 10))[0]
        x = sm.add_constant(a_values[mask])
        y = m_values[mask]
        self.ma_model = sm.RLM(y, x).fit().params

    def normalize(self):
        for chrom in self.peaks1.peaks:
            for peak in self.peaks1.peaks[chrom]:
                peak.normalize(self.ma_model)
        for chrom in self.peaks2.peaks:
            for peak in self.peaks2.peaks[chrom]:
                peak.normalize(self.ma_model)
        for chrom in self.common_peaks.peaks:
            for peak in self.common_peaks.peaks[chrom]:
                peak.normalize(self.ma_model)
