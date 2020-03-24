"""
manorm.io
---------

This module contains the input/output related functions.
"""

import os
from collections import defaultdict
from math import log10


def mk_dir(root_dir):
    if not os.path.isdir(root_dir):
        os.makedirs(root_dir)
    output_dirs = {'figures': os.path.join(root_dir, 'output_figures'),
                   'filters': os.path.join(root_dir, 'output_filters'),
                   'tracks': os.path.join(root_dir, 'output_tracks')}
    for sub_dir in output_dirs.values():
        if not os.path.isdir(sub_dir):
            os.makedirs(sub_dir)


def _get_unique_and_merged_peaks(peaks1, peaks2, peaks_merged):
    peaks = []
    peak_groups = []
    for chrom in peaks1.chroms:
        for peak in peaks1.fetch(chrom):
            if not peak.iscommon:
                peaks.append(peak)
                peak_groups.append(peaks1.name + '_unique')
    for chrom in peaks_merged.chroms:
        for peak in peaks_merged.fetch(chrom):
            peaks.append(peak)
            peak_groups.append('merged_common')
    for chrom in peaks2.chroms:
        for peak in peaks2.fetch(chrom):
            if not peak.iscommon:
                peaks.append(peak)
                peak_groups.append(peaks2.name + '_unique')
    return peaks, peak_groups


def write_original_peaks(root_dir, peaks1, peaks2):
    sample_names = [peaks1.name, peaks2.name]
    peaks = [peaks1, peaks2]

    header = f"chr\tstart\tend\tsummit\tM_value\tA_value\tP_value\t" \
             f"Peak_Group\tnormalized_read_density_in_{peaks1.name}\t" \
             f"normalized_read_density_in_{peaks2.name}\n"

    for temp_name, temp_peaks in zip(sample_names, peaks):
        temp_file = os.path.join(root_dir, temp_name + '_MAvalues.xls')
        with open(temp_file, 'w') as fout:
            fout.write(header)
            for chrom in temp_peaks.chroms:
                for peak in temp_peaks.fetch(chrom):
                    if peak.iscommon:
                        peak_group = temp_name + "_common"
                    else:
                        peak_group = temp_name + "_unique"
                    fout.write(
                        f"{peak.chrom}\t{peak.start + 1}\t{peak.end}\t"
                        f"{peak.summit + 1}\t{peak.m_normed:.5f}\t"
                        f"{peak.a_normed:.5f}\t{peak.p_value}\t{peak_group}\t"
                        f"{peak.read_density1_normed:.5f}\t"
                        f"{peak.read_density2_normed:.5f}\n")


def write_all_peaks(root_dir, peaks1, peaks2, peaks_merged):
    peaks, peak_groups = _get_unique_and_merged_peaks(peaks1, peaks2,
                                                      peaks_merged)
    header = f"chr\tstart\tend\tsummit\tM_value\tA_value\tP_value\t" \
             f"Peak_Group\tnormalized_read_density_in_{peaks1.name}\t" \
             f"normalized_read_density_in_{peaks2.name}\n"
    path = os.path.join(
        root_dir, peaks1.name + '_vs_' + peaks2.name + '_all_MAvalues.xls')
    with open(path, 'w') as fout:
        fout.write(header)
        for peak, peak_group in zip(peaks, peak_groups):
            fout.write(
                f"{peak.chrom}\t{peak.start + 1}\t{peak.end}\t"
                f"{peak.summit + 1}\t{peak.m_normed:.5f}\t"
                f"{peak.a_normed:.5f}\t{peak.p_value}\t{peak_group}\t"
                f"{peak.read_density1_normed:.5f}\t"
                f"{peak.read_density2_normed:.5f}\n")


def write_wiggle_track(root_dir, peaks1, peaks2, peaks_merged):
    peaks, _ = _get_unique_and_merged_peaks(peaks1, peaks2, peaks_merged)
    tracks = defaultdict(list)
    for peak in peaks:
        tracks[peak.chrom].append((peak.summit + 1, peak.m_normed,
                                   peak.a_normed, peak.p_value))
    for chrom in tracks:
        tracks[chrom].sort(key=lambda x: x[0])

    # write tracks for M values, A values and P values
    output_prefix = peaks1.name + '_vs_' + peaks2.name
    path_m = os.path.join(
        root_dir, 'output_tracks', output_prefix + '_M_values.wig')
    path_a = os.path.join(
        root_dir, 'output_tracks', output_prefix + '_A_values.wig')
    path_p = os.path.join(
        root_dir, 'output_tracks', output_prefix + '_P_values.wig')
    with open(path_m, 'w') as fout_m, open(path_a, 'w') as fout_a, \
            open(path_p, 'w') as fout_p:
        fout_m.write(
            f"track type=wiggle_0 name={output_prefix}_M_value "
            f"visibility=full autoScale=on color=255,0,0 yLineMark=0 "
            f"yLineOnOff=on priority=10\n")
        fout_a.write(
            f"track type=wiggle_0 name={output_prefix}_A_value "
            f"visibility=full autoScale=on color=255,0,0 yLineMark=0 "
            f"yLineOnOff=on priority=10\n")
        fout_p.write(
            f"track type=wiggle_0 name={output_prefix}_-log10(P_value) "
            f"visibility=full autoScale=on color=255,0,0 yLineMark=0 "
            f"yLineOnOff=on priority=10\n")
        for chrom in tracks:
            if len(tracks[chrom]) == 0:
                continue
            fout_m.write(f"variableStep chrom={chrom} span=100\n")
            fout_a.write(f"variableStep chrom={chrom} span=100\n")
            fout_p.write(f"variableStep chrom={chrom} span=100\n")
            for summit, m_value, a_value, p_value in tracks[chrom]:
                fout_m.write(f"{summit}\t{m_value:.5f}\n")
                fout_a.write(f"{summit}\t{a_value:.5f}\n")
                fout_p.write(f"{summit}\t{-log10(p_value)}\n")


def write_biased_peaks(root_dir, peaks1, peaks2, peaks_merged, m_cutoff,
                       p_cutoff):
    m_cutoff = abs(m_cutoff)
    peaks, peak_groups = _get_unique_and_merged_peaks(peaks1, peaks2,
                                                      peaks_merged)
    output_prefix = peaks1.name + '_vs_' + peaks2.name
    path_biased1 = os.path.join(root_dir, 'output_filters', output_prefix
                                + f"_M_above_{m_cutoff}_biased_peaks.bed")
    path_biased2 = os.path.join(root_dir, 'output_filters', output_prefix
                                + f"_M_below_-{m_cutoff}_biased_peaks.bed")
    path_unbiased = os.path.join(root_dir, 'output_filters',
                                 output_prefix + '_unbiased_peaks.bed')
    num_biased1, num_biased2, num_unbiased = 0, 0, 0
    with open(path_biased1, 'w') as fout_biased1, \
            open(path_biased2, 'w') as fout_biased2, \
            open(path_unbiased, 'w') as fout_unbiased:
        for peak, peak_group in zip(peaks, peak_groups):
            line = f"{peak.chrom}\t{peak.start}\t{peak.end}\t{peak_group}\t" \
                   f"{peak.m_normed:.5f}\n"
            if abs(peak.m_normed) < m_cutoff:
                num_unbiased += 1
                fout_unbiased.write(line)
            elif peak.p_value <= p_cutoff:
                if peak.m_normed >= m_cutoff:
                    num_biased1 += 1
                    fout_biased1.write(line)
                elif peak.m_normed <= -m_cutoff:
                    num_biased2 += 1
                    fout_biased2.write(line)
    return num_biased1, num_biased2, num_unbiased
