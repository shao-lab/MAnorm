"""Core module for MAnorm."""

from collections import defaultdict
import numpy as np
import statsmodels.api as sm
from manorm.lib.peak import Peaks, overlap_on_single_chr, generate_random_peaks, merge_peaks
from manorm.lib.read import Reads


class MAnorm(object):
    def __init__(self, name1, name2, root_dir, output_prefix):
        self.name1 = name1
        self.name2 = name2
        self.peaks1 = Peaks()
        self.peaks2 = Peaks()
        self.common_peaks = Peaks()
        self.common_dist = defaultdict(list)
        self.reads1 = Reads()
        self.reads2 = Reads()
        self.ma_model = None
        self.root_dir = root_dir
        if self.name1 and self.name2:
            self.output_prefix = self.name1 + "_vs_" + self.name2
        else:
            self.output_prefix = output_prefix

    def load_peaks(self, peaks_file1, peaks_file2):
        self.peaks1 = Peaks(peaks_file1, self.name1)
        self.peaks2 = Peaks(peaks_file2, self.name2)

    def load_reads(self, reads_file1, reads_file2, shift_size1, shift_size2):
        self.reads1 = Reads(reads_file1, self.name1, shift_size1)
        self.reads2 = Reads(reads_file2, self.name2, shift_size2)

    def classify_peaks_by_overlap(self):
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
