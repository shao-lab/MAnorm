# -*- coding: utf-8 -*-

"""
manorm.io
~~~~~~~~~

This module contains the input/output related functions.
"""

from __future__ import absolute_import

import os
from collections import defaultdict
from math import log10

from manorm.compat import zip


def mk_dir(root_dir):
    if not os.path.isdir(root_dir):
        os.makedirs(root_dir)
    output_dirs = {'figures': os.path.join(root_dir, 'output_figures'),
                   'filters': os.path.join(root_dir, 'output_filters'),
                   'tracks': os.path.join(root_dir, 'output_tracks')}
    for key in output_dirs:
        if not os.path.isdir(output_dirs[key]):
            os.makedirs(output_dirs[key])


def write_original_peaks(root_dir, peaks1, peaks2):
    sample_names = [peaks1.name, peaks2.name]
    peaks = [peaks1, peaks2]

    header = "chr\tstart\tend\tsummit\tM_value\tA_value\tP_value\tPeak_Group\tnormalized_read_density_in_{}\t" \
             "normalized_read_density_in_{}\n".format(peaks1.name, peaks2.name)
    formatter = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n"

    for temp_name, temp_peaks in zip(sample_names, peaks):
        temp_file = os.path.join(root_dir, temp_name + '_MAvalues.xls')
        with open(temp_file, 'w') as fout:
            fout.write(header)
            for chrom in temp_peaks.chroms:
                for peak in temp_peaks.fetch(chrom):
                    fout.write(formatter.format(peak.chrom, peak.start + 1, peak.end, peak.summit - peak.start,
                                                peak.m_value_normed, peak.a_value_normed, peak.p_value,
                                                temp_name + '_' + peak.type, peak.read_density1_normed,
                                                peak.read_density2_normed))


def write_all_peaks(root_dir, peaks1, peaks2, peaks_merged):
    file_name = os.path.join(root_dir, peaks1.name + '_vs_' + peaks2.name + '_all_MAvalues.xls')
    peaks = []
    peak_types = []
    for chrom in peaks1.chroms:
        for peak in peaks1.fetch(chrom):
            if peak.type == 'unique':
                peaks.append(peak)
                peak_types.append(peaks1.name + '_unique')
    for chrom in peaks_merged.chroms:
        for peak in peaks_merged.fetch(chrom):
            peaks.append(peak)
            peak_types.append('merged_common')

    for chrom in peaks2.chroms:
        for peak in peaks2.fetch(chrom):
            if peak.type == 'unique':
                peaks.append(peak)
                peak_types.append(peaks2.name + '_unique')

    header = "chr\tstart\tend\tsummit\tM_value\tA_value\tP_value\tPeak_Group\tnormalized_read_density_in_{}\t" \
             "normalized_read_density_in_{}\n".format(peaks1.name, peaks2.name)
    formatter = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n"
    with open(file_name, 'w') as fout:
        fout.write(header)
        for idx, peak in enumerate(peaks):
            fout.write(formatter.format(peak.chrom, peak.start + 1, peak.end, peak.summit - peak.start,
                                        peak.m_value_normed, peak.a_value_normed, peak.p_value, peak_types[idx],
                                        peak.read_density1_normed, peak.read_density2_normed))


def write_wiggle_track(root_dir, peaks1, peaks2, peaks_merged):
    output_prefix = peaks1.name + '_vs_' + peaks2.name
    tracks = defaultdict(list)
    for chrom in peaks1.chroms:
        for peak in peaks1.fetch(chrom):
            if peak.type == 'unique':
                tracks[chrom].append((peak.summit + 1, peak.m_value_normed, peak.a_value_normed, peak.p_value))
    for chrom in peaks_merged.chroms:
        for peak in peaks_merged.fetch(chrom):
            tracks[chrom].append((peak.summit + 1, peak.m_value_normed, peak.a_value_normed, peak.p_value))
    for chrom in peaks2.chroms:
        for peak in peaks2.fetch(chrom):
            if peak.type == 'unique':
                tracks[chrom].append((peak.summit + 1, peak.m_value_normed, peak.a_value_normed, peak.p_value))
    for chrom in tracks:
        tracks[chrom].sort(key=lambda x: x[0])

    file_name = os.path.join(root_dir, 'output_tracks', output_prefix + '_M_values.wig')
    with open(file_name, 'w') as fout:
        fout.write("track type=wiggle_0 name={}_M_value visibility=full autoScale=on color=255,0,0 yLineMark=0 "
                   "yLineOnOff=on priority=10\n".format(output_prefix))
        for chrom in tracks:
            if len(tracks[chrom]) == 0:
                continue
            fout.write("variableStep chrom={} span=100\n".format(chrom))
            for summit, m_value, _, _ in tracks[chrom]:
                fout.write("{}\t{}\n".format(summit, m_value))

    file_name = os.path.join(root_dir, 'output_tracks', output_prefix + '_A_values.wig')
    with open(file_name, 'w') as fout:
        fout.write("track type=wiggle_0 name={}_A_value visibility=full autoScale=on color=255,0,0 yLineMark=0 "
                   "yLineOnOff=on priority=10\n".format(output_prefix))
        for chrom in tracks:
            if len(tracks[chrom]) == 0:
                continue
            fout.write("variableStep chrom={} span=100\n".format(chrom))
            for summit, _, a_value, _ in tracks[chrom]:
                fout.write("{}\t{}\n".format(summit, a_value))

    file_name = os.path.join(root_dir, 'output_tracks', output_prefix + '_P_values.wig')
    with open(file_name, 'w') as fout:
        fout.write("track type=wiggle_0 name={}_-log10(P_value) visibility=full autoScale=on color=255,0,0 yLineMark=0 "
                   "yLineOnOff=on priority=10\n".format(output_prefix))
        for chrom in tracks:
            if len(tracks[chrom]) == 0:
                continue
            fout.write("variableStep chrom={} span=100\n".format(chrom))
            for summit, _, _, p_value in tracks[chrom]:
                fout.write("{}\t{}\n".format(summit, -log10(p_value)))


def write_biased_peaks(root_dir, peaks1, peaks2, peaks_merged, m_cutoff, p_cutoff):
    m_cutoff = abs(m_cutoff)
    peaks = defaultdict(list)
    peak_types = defaultdict(list)
    for chrom in peaks1.chroms:
        for peak in peaks1.fetch(chrom):
            if peak.type == 'unique':
                peaks[chrom].append(peak)
                peak_types[chrom].append(peaks1.name + '_unique')
    for chrom in peaks_merged.chroms:
        for peak in peaks_merged.fetch(chrom):
            peaks[chrom].append(peak)
            peak_types[chrom].append('merged_common')
    for chrom in peaks2.chroms:
        for peak in peaks2.fetch(chrom):
            if peak.type == 'unique':
                peaks[chrom].append(peak)
                peak_types[chrom].append(peaks2.name + '_unique')

    output_prefix = peaks1.name + '_vs_' + peaks2.name
    num_biased1, num_biased2, num_unbiased = 0, 0, 0
    file_name1 = os.path.join(root_dir, 'output_filters',
                              output_prefix + "_M_above_{}_biased_peaks.bed".format(m_cutoff))
    file_name2 = os.path.join(root_dir, 'output_filters',
                              output_prefix + "_M_below_-{}_biased_peaks.bed".format(m_cutoff))
    file_name = os.path.join(root_dir, 'output_filters', output_prefix + '_unbiased_peaks.bed')
    with open(file_name1, 'w') as fout1, open(file_name2, 'w') as fout2, open(file_name, 'w') as fout3:
        for chrom in peaks:
            for idx, peak in enumerate(peaks[chrom]):
                if abs(peak.m_value_normed) < m_cutoff:
                    num_unbiased += 1
                    fout3.write("{}\t{}\t{}\t{}\t{}\n".format(peak.chrom, peak.start, peak.end,
                                                              peak_types[chrom][idx], peak.m_value_normed))
                elif peak.p_value <= p_cutoff:
                    if peak.m_value_normed >= m_cutoff:
                        num_biased1 += 1
                        fout1.write("{}\t{}\t{}\t{}\t{}\n".format(peak.chrom, peak.start, peak.end,
                                                                  peak_types[chrom][idx], peak.m_value_normed))
                    elif peak.m_value_normed <= -m_cutoff:
                        num_biased2 += 1
                        fout2.write("{}\t{}\t{}\t{}\t{}\n".format(peak.chrom, peak.start, peak.end,
                                                                  peak_types[chrom][idx], peak.m_value_normed))
    return num_biased1, num_biased2, num_unbiased
