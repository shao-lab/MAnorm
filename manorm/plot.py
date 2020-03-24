"""
manorm.plot
-----------

This module contains the plot functions.
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
            if not peak.iscommon:
                peaks1_unique.append(peak)
    peaks2_unique = []
    for chrom in peaks2.chroms:
        for peak in peaks2.fetch(chrom):
            if not peak.iscommon:
                peaks2_unique.append(peak)
    merged_common_peaks = []
    for chrom in peaks_merged.chroms:
        for peak in peaks_merged.fetch(chrom):
            merged_common_peaks.append(peak)

    output_prefix = peaks1.name + '_vs_' + peaks2.name

    peaks1_name = peaks1.name + '_unique'
    peaks2_name = peaks2.name + '_unique'
    merged_peaks_name = 'merged_common_peaks'
    peaks_names = [peaks1_name, peaks2_name, merged_peaks_name]
    colors = ["#E53A40", "#30A9DE", "#566270"]

    # plot the relationship of read densities
    fig, ax = plt.subplots(figsize=(4, 4))
    reads_density1, reads_density2 = [], []
    for peak in merged_common_peaks:
        reads_density1.append(peak.read_density1)
        reads_density2.append(peak.read_density2)
    x = np.log2(reads_density1)
    y = np.log2(reads_density2)
    x_max = max(x)
    x_min = min(x)
    ax.scatter(x, y, s=1, c="#566270", label=merged_peaks_name, alpha=0.8)
    rx = np.arange(x_min, x_max, 0.01)
    ry = (2 - ma_params[1]) * rx / (2 + ma_params[1]) - 2 * ma_params[0] / (
            2 + ma_params[1])
    ax.plot(rx, ry, "-", color="#1E2022")
    ax.legend(loc='upper left', fontsize=6, handletextpad=0, markerscale=2,
              frameon=False)
    ax.tick_params(labelsize=8, pad=2)
    ax.set_xlabel(f"$log_2$ read density in {peaks1.name}", fontsize=8)
    ax.set_ylabel(f"$log_2$ read density in {peaks2.name}", fontsize=8)
    ax.set_title("M-A model fitted on common peaks", fontsize=10)
    fig.subplots_adjust(left=0.15, right=0.9, bottom=0.15, top=0.9)
    plt.savefig(os.path.join(root_dir, 'output_figures', output_prefix +
                             '_read_density_on_common_peaks.pdf'))

    # plot the MA plot before normalization
    fig, ax = plt.subplots(figsize=(4, 3))
    a_max = 0
    a_min = 999999999
    for idx, peaks in enumerate(
            [peaks1_unique, peaks2_unique, merged_common_peaks]):
        m_values = [peak.m_raw for peak in peaks]
        a_values = [peak.a_raw for peak in peaks]
        a_max = max(max(a_values), a_max)
        a_min = min(min(a_values), a_min)
        plt.scatter(a_values, m_values, s=1, c=colors[idx],
                    label=peaks_names[idx], alpha=0.8)
    ax.axhline(y=0, ls='--', color='lightgrey')
    x = np.arange(a_min, a_max, 0.01)
    y = ma_params[1] * x + ma_params[0]
    ax.plot(x, y, "-", color="#1E2022")
    ymin, ymax = ax.get_ylim()
    ylim = max(abs(ymin), abs(ymax))
    ax.set_ylim(ymin=-ylim, ymax=ylim)
    ax.tick_params(labelsize=8, pad=2)
    ax.set_xlabel("$A$ value", fontsize=8, labelpad=3)
    ax.set_ylabel("$M$ value", fontsize=8, labelpad=3)
    ax.set_title("M-A plot before normalization", fontsize=10)
    ax.legend(loc='upper right', fontsize=6, scatteryoffsets=[0.5],
              handletextpad=0, markerscale=2, frameon=False)
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.9)
    fig.savefig(os.path.join(root_dir, 'output_figures', output_prefix +
                             '_MA_plot_before_normalization.pdf'))

    # plot the MA plot after normalization
    fig, ax = plt.subplots(figsize=(4, 3))
    for idx, peaks in enumerate(
            [peaks1_unique, peaks2_unique, merged_common_peaks]):
        m_values = [peak.m_normed for peak in peaks]
        a_values = [peak.a_normed for peak in peaks]
        plt.scatter(a_values, m_values, s=1, c=colors[idx],
                    label=peaks_names[idx], alpha=0.8)
    ax.axhline(y=0, ls='--', color='lightgrey')
    ymin, ymax = ax.get_ylim()
    ylim = max(abs(ymin), abs(ymax))
    ax.set_ylim(ymin=-ylim, ymax=ylim)
    ax.tick_params(labelsize=8, pad=2)
    ax.set_xlabel("Normalized $A$ value", fontsize=8, labelpad=3)
    ax.set_ylabel("Normalized $M$ value", fontsize=8, labelpad=3)
    ax.set_title("M-A plot after normalization", fontsize=10)
    ax.legend(loc='upper right', fontsize=6, scatteryoffsets=[0.5],
              handletextpad=0, markerscale=2, frameon=False)
    fig.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.9)
    fig.savefig(os.path.join(root_dir, 'output_figures', output_prefix +
                             '_MA_plot_after_normalization.pdf'))

    # plot the MA plot after normalization colored by P value
    fig, ax = plt.subplots(figsize=(4, 3))
    m_values = []
    a_values = []
    p_values = []
    for idx, peaks in enumerate(
            [peaks1_unique, peaks2_unique, merged_common_peaks]):
        for peak in peaks:
            m_values.append(peak.m_normed)
            a_values.append(peak.a_normed)
            p_values.append(peak.p_value)
    colors = -np.log10(p_values)
    for i, c in enumerate(colors):
        if c > 50:
            colors[i] = 50
    scatter = ax.scatter(a_values, m_values, s=1, c=colors, cmap="coolwarm")
    ax.axhline(y=0, ls='--', color='lightgrey')
    ymin, ymax = ax.get_ylim()
    ylim = max(abs(ymin), abs(ymax))
    ax.set_ylim(ymin=-ylim, ymax=ylim)
    ax.tick_params(labelsize=8, pad=2)
    ax.set_xlabel("Normalized $A$ value", fontsize=8, labelpad=3)
    ax.set_ylabel("Normalized $M$ value", fontsize=8, labelpad=3)
    ax.set_title("M-A plot colored by $P$-value", fontsize=10)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.ax.tick_params(labelsize=6, pad=2)
    cbar.set_label("$-log_{10}$($P$-value)", fontsize=7)
    fig.subplots_adjust(left=0.15, right=0.98, bottom=0.15, top=0.9)
    fig.savefig(os.path.join(root_dir, 'output_figures', output_prefix +
                             '_MA_plot_with_P_value.pdf'))
