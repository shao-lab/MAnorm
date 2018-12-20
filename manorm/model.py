# -*- coding: utf-8 -*-

"""
manorm.model
~~~~~~~~~~~~

This module contains classes and functions designed for MAnorm normalization model.
"""

from __future__ import absolute_import

import numpy as np
import statsmodels.api as sm

from manorm.peak.utils import classify_peaks_by_overlap, merge_common_peaks


class MAmodel(object):
    def __init__(self, peaks1, peaks2, reads1, reads2):
        self.peaks1 = peaks1
        self.peaks2 = peaks2
        self.peaks_merged = None
        self.reads1 = reads1
        self.reads2 = reads2
        self.ma_params = None

    def process_peaks(self):
        """Preprocess peaks."""
        self.peaks1, self.peaks2 = classify_peaks_by_overlap(self.peaks1, self.peaks2)
        self.peaks_merged = merge_common_peaks(self.peaks1, self.peaks2)

    def cal_ma_value(self, window_size=2000):
        """Calculate m values and a values of peaks."""
        for chrom in self.peaks1.chroms:
            for peak in self.peaks1.fetch(chrom):
                peak.count_reads(self.reads1, self.reads2, window_size)
                peak.cal_ma_value()
        for chrom in self.peaks2.chroms:
            for peak in self.peaks2.fetch(chrom):
                peak.count_reads(self.reads1, self.reads2, window_size)
                peak.cal_ma_value()
        for chrom in self.peaks_merged.chroms:
            for peak in self.peaks_merged.fetch(chrom):
                peak.count_reads(self.reads1, self.reads2, window_size)
                peak.cal_ma_value()

    def fit_model(self, window_size=2000, summit_dis=500):
        """Fit M-A normalization model."""
        self.cal_ma_value(window_size=window_size)
        m_values = []
        a_values = []
        for chrom in self.peaks_merged.chroms:
            for idx, peak in enumerate(self.peaks_merged.fetch(chrom)):
                if peak.summit_dis <= summit_dis:
                    m_values.append(peak.m_value)
                    a_values.append(peak.a_value)
        m_values = np.array(m_values)
        a_values = np.array(a_values)
        mask = np.where((m_values >= -10) & (m_values <= 10))[0]
        x = sm.add_constant(a_values[mask])
        y = m_values[mask]
        self.ma_params = sm.RLM(y, x).fit().params

    def normalize(self):
        """Normalize all peaks."""
        for chrom in self.peaks1.chroms:
            for peak in self.peaks1.fetch(chrom):
                peak.normalize(self.ma_params)
        for chrom in self.peaks2.chroms:
            for peak in self.peaks2.fetch(chrom):
                peak.normalize(self.ma_params)
        for chrom in self.peaks_merged.chroms:
            for peak in self.peaks_merged.fetch(chrom):
                peak.normalize(self.ma_params)
