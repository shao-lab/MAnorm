# -*- coding: utf-8 -*-

"""
manorm.plot
~~~~~~~~~~~

This module contains plot functions.
"""

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def plt_figures(root_dir, peaks1, peaks2, peaks_merged, ma_params):
    peaks1_unique = []
    for chrom in peaks1.chroms:
        for peak in peaks1.fetch(chrom):
            if peak.type == 'unique':
                peaks1_unique.append(peak)

    peaks2_unique = []
    for chrom in peaks2.chroms:
        for peak in peaks2.fetch(chrom):
            if peak.type == 'unique':
                peaks2_unique.append(peak)

    merged_common_peaks = []
    for chrom in peaks_merged.chroms:
        for peak in peaks_merged.fetch(chrom):
            merged_common_peaks.append(peak)

    peaks1_name = peaks1.name + '_unique'
    peaks2_name = peaks2.name + '_unique'
    merged_peaks_name = 'merged_common_peaks'

    peaks_names = [peaks1_name, peaks2_name, merged_peaks_name]
    colors = ["#E53A40", "#30A9DE", "#566270"]

    output_prefix = peaks1.name + '_vs_' + peaks2.name
    fig, ax = plt.subplots(figsize=(8, 6))
    reads_density1, reads_density2 = [], []
    for peak in merged_common_peaks:
        reads_density1.append(peak.read_density1)
        reads_density2.append(peak.read_density2)
    log_read_density_max = max(np.log2(reads_density1))
    log_read_density_min = min(np.log2(reads_density1))
    ax.scatter(np.log2(reads_density1), np.log2(reads_density2), s=6, c="#566270", label=merged_peaks_name, alpha=0.5)
    rx = np.arange(log_read_density_min, log_read_density_max, 0.01)
    ry = (2 - ma_params[1]) * rx / (2 + ma_params[1]) - 2 * ma_params[0] / (2 + ma_params[1])
    ax.plot(rx, ry, "-", color="#1E2022")
    ax.legend(loc=1, fontsize=8)
    ax.set_xlabel("log2 read density in {}".format(peaks1.name))
    ax.set_ylabel("log2 read density in {}".format(peaks2.name))
    ax.set_title("M-A model fitted via common peaks", fontsize=18)
    plt.savefig(os.path.join(root_dir, 'output_figures', output_prefix + '_read_density_on_common_peaks.png'),
                dpi=300)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.grid(axis="y", linestyle="--")
    a_max = 0
    a_min = 999999999
    for idx, peaks in enumerate([peaks1_unique, peaks2_unique, merged_common_peaks]):
        m_values = [peak.m_value for peak in peaks]
        a_values = [peak.a_value for peak in peaks]
        a_max = max(max(a_values), a_max)
        a_min = min(min(a_values), a_min)
        plt.scatter(a_values, m_values, s=6, c=colors[idx], label=peaks_names[idx], alpha=0.5)
    ax.legend(loc=1, fontsize=8)

    x = np.arange(a_min, a_max, 0.01)
    y = ma_params[1] * x + ma_params[0]
    ax.plot(x, y, "-", color="#1E2022")
    ax.set_xlabel("A value", fontsize=16)
    ax.set_ylabel("M value", fontsize=16)
    ax.set_title("M-A plot before normalization", fontsize=18)
    fig.savefig(os.path.join(root_dir, 'output_figures', output_prefix + '_MA_plot_before_normalization.png'),
                dpi=300)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.grid(axis="y", linestyle="--")
    for idx, peaks in enumerate([peaks1_unique, peaks2_unique, merged_common_peaks]):
        m_values = [peak.m_value_normed for peak in peaks]
        a_values = [peak.a_value_normed for peak in peaks]
        plt.scatter(a_values, m_values, s=6, c=colors[idx], label=peaks_names[idx], alpha=0.5)
    ax.legend(loc=1, fontsize=8)
    ax.set_xlabel("A value", fontsize=16)
    ax.set_ylabel("M value", fontsize=16)
    ax.set_title("M-A plot after normalization", fontsize=18)
    fig.savefig(os.path.join(root_dir, 'output_figures', output_prefix + '_MA_plot_after_normalization.png'),
                dpi=300)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.grid(axis="y", linestyle="--")
    m_values = []
    a_values = []
    p_values = []
    for idx, peaks in enumerate([peaks1_unique, peaks2_unique, merged_common_peaks]):
        for peak in peaks:
            m_values.append(peak.m_value_normed)
            a_values.append(peak.a_value_normed)
            p_values.append(peak.p_value)
    colors = -np.log10(p_values)
    for i, c in enumerate(colors):
        if c > 50:
            colors[i] = 50
    scatter = ax.scatter(a_values, m_values, s=10, c=colors, cmap="coolwarm")
    fig.colorbar(scatter, ax=ax)
    ax.set_xlabel("A value", fontsize=16)
    ax.set_ylabel("M value", fontsize=16)
    ax.set_title("-log10(P-value)", fontsize=18)
    fig.savefig(os.path.join(root_dir, 'output_figures', output_prefix + '_MA_plot_with_P_value.png'), dpi=300)
