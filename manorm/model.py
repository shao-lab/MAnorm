"""
manorm.model
------------

This module contains the core MAnorm model.
"""

import numpy as np
from sklearn.linear_model import HuberRegressor

from manorm.exceptions import ProcessNotReadyError
from manorm.region.utils import classify_peaks_by_overlap, merge_common_peaks


class MAmodel(object):
    def __init__(self, peaks1, peaks2, reads1, reads2):
        self.peaks1 = peaks1
        self.peaks2 = peaks2
        self.peaks_merged = None
        self.reads1 = reads1
        self.reads2 = reads2
        self.ma_params = None
        self.processed = False
        self.fitted = False
        self.normalized = False

    def process_peaks(self):
        """Preprocess peaks."""
        self.peaks1, self.peaks2 = classify_peaks_by_overlap(self.peaks1,
                                                             self.peaks2)
        self.peaks_merged = merge_common_peaks(self.peaks1, self.peaks2)
        self.processed = True

    def _count_reads(self, window_size=2000):
        """Calculate m values and a values of peaks."""
        for chrom in self.peaks1.chroms:
            for peak in self.peaks1.fetch(chrom):
                peak.count_reads(self.reads1, self.reads2, window_size)
        for chrom in self.peaks2.chroms:
            for peak in self.peaks2.fetch(chrom):
                peak.count_reads(self.reads1, self.reads2, window_size)
        for chrom in self.peaks_merged.chroms:
            for peak in self.peaks_merged.fetch(chrom):
                peak.count_reads(self.reads1, self.reads2, window_size)

    def fit_model(self, window_size=2000, summit_dis_cutoff=500):
        """Fit M-A normalization model."""
        if not self.processed:
            raise ProcessNotReadyError("fit the M-A model", 'process peaks')
        self._count_reads(window_size=window_size)
        m_values = []
        a_values = []
        for chrom in self.peaks_merged.chroms:
            for peak in self.peaks_merged.fetch(chrom):
                if peak.summit_dis <= summit_dis_cutoff:
                    m_values.append(peak.m_raw)
                    a_values.append(peak.a_raw)
        m_values = np.array(m_values)
        a_values = np.array(a_values)
        mask = abs(m_values) <= 10
        huber = HuberRegressor()
        huber.fit(a_values[mask].reshape(-1, 1), m_values[mask])
        self.ma_params = [huber.intercept_, huber.coef_[0]]
        self.fitted = True

    def normalize(self):
        """Normalize all peaks."""
        if not self.fitted:
            raise ProcessNotReadyError("normalize peaks", "fit the M-A model")
        intercept = self.ma_params[0]
        slope = self.ma_params[1]
        for chrom in self.peaks1.chroms:
            for peak in self.peaks1.fetch(chrom):
                peak.normalize(slope=slope, intercept=intercept)
        for chrom in self.peaks2.chroms:
            for peak in self.peaks2.fetch(chrom):
                peak.normalize(slope=slope, intercept=intercept)
        for chrom in self.peaks_merged.chroms:
            for peak in self.peaks_merged.fetch(chrom):
                peak.normalize(slope=slope, intercept=intercept)
        self.normalized = True
