"""IO module of MAnorm."""

import os
from math import log10
from collections import defaultdict

__all__ = ["mk_dir",  "output_original_peaks","output_all_peaks", "output_wiggle_track",
           "output_unbiased_peaks", "output_biased_peaks"]


def mk_dir(root_dir):
    if not os.path.isdir(root_dir):
        os.makedirs(root_dir)

    output_dirs = {"figures": os.path.join(root_dir, "output_figures"),
                   "filters": os.path.join(root_dir, "output_filters"),
                   "tracks": os.path.join(root_dir, "output_tracks")}

    for key in output_dirs:
        if not os.path.isdir(output_dirs[key]):
            os.makedirs(output_dirs[key])


def output_original_peaks(ma):
    file_name1 = os.path.join(ma.root_dir, ma.peaks1.name + "_MAvalues.xls")
    file_name2 = os.path.join(ma.root_dir, ma.peaks2.name + "_MAvalues.xls")
    file_names = [file_name1, file_name2]
    peaks = [ma.peaks1.peaks, ma.peaks2.peaks]

    header = "chr\tstart\tend\tsummit\tM_value\tA_value\tP_value\tPeak_Group\tnormalized_read_density_in_{}\t" \
             "normalized_read_density_in_{}\n".format(ma.reads1.name, ma.reads2.name)
    formatter = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n"

    for temp_file, temp_peaks in zip(file_names, peaks):
        with open(temp_file, "w") as fout:
            fout.write(header)
            for chrom in temp_peaks:
                for peak in temp_peaks[chrom]:
                    fout.write(formatter.format(peak.chrom, peak.start + 1, peak.end, peak.summit - peak.start,
                                                peak.normed_m_value, peak.normed_a_value, peak.p_value,
                                                peak.type, peak.normed_read_density1, peak.normed_read_density2))


def output_all_peaks(ma):
    file_name = os.path.join(ma.root_dir, ma.output_prefix + "_all_MAvalues.xls")
    peaks = []
    for chrom in ma.peaks1.peaks:
        for peak in ma.peaks1.peaks[chrom]:
            if peak.type.endswith("unique"):
                peaks.append(peak)
    for chrom in ma.common_peaks.peaks:
        for peak in ma.common_peaks.peaks[chrom]:
            peaks.append(peak)
    for chrom in ma.peaks2.peaks:
        for peak in ma.peaks2.peaks[chrom]:
            if peak.type.endswith("unique"):
                peaks.append(peak)

    header = "chr\tstart\tend\tsummit\tM_value\tA_value\tP_value\tPeak_Group\tnormalized_read_density_in_{}\t" \
             "normalized_read_density_in_{}\n".format(ma.reads1.name, ma.reads2.name)
    formatter = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n"
    with open(file_name, "w") as fout:
        fout.write(header)
        for peak in peaks:
            fout.write(formatter.format(peak.chrom, peak.start + 1, peak.end, peak.summit - peak.start,
                                        peak.normed_m_value, peak.normed_a_value, peak.p_value,
                                        peak.type, peak.normed_read_density1, peak.normed_read_density2))


def output_wiggle_track(ma):
    tracks = defaultdict(list)
    for chrom in ma.peaks1.peaks:
        for peak in ma.peaks1.peaks[chrom]:
            if peak.type.endswith("unique"):
                tracks[chrom].append((peak.summit + 1, peak.normed_m_value, peak.normed_a_value, peak.p_value))
    for chrom in ma.common_peaks.peaks:
        for peak in ma.common_peaks.peaks[chrom]:
            tracks[chrom].append((peak.summit + 1, peak.normed_m_value, peak.normed_a_value, peak.p_value))
    for chrom in ma.peaks2.peaks:
        for peak in ma.peaks2.peaks[chrom]:
            if peak.type.endswith("unique"):
                tracks[chrom].append((peak.summit + 1, peak.normed_m_value, peak.normed_a_value, peak.p_value))
    for chrom in tracks:
        tracks[chrom].sort(key=lambda x: x[0])

    file_name = os.path.join(ma.root_dir, "output_tracks", ma.output_prefix + "_M_values.wig")
    with open(file_name, "w") as fout:
        fout.write("track type=wiggle_0 name={}_M_value visibility=full autoScale=on color=255,0,0 yLineMark=0 "
                   "yLineOnOff=on priority=10\n".format(ma.output_prefix))
        for chrom in tracks:
            if len(tracks[chrom]) == 0:
                continue
            fout.write("variableStep chrom={} span=100\n".format(chrom))
            for summit, m_value, _, _ in tracks[chrom]:
                fout.write("{}\t{}\n".format(summit, m_value))

    file_name = os.path.join(ma.root_dir, "output_tracks", ma.output_prefix + "_A_values.wig")
    with open(file_name, "w") as fout:
        fout.write("track type=wiggle_0 name={}_A_value visibility=full autoScale=on color=255,0,0 yLineMark=0 "
                   "yLineOnOff=on priority=10\n".format(ma.output_prefix))
        for chrom in tracks:
            if len(tracks[chrom]) == 0:
                continue
            fout.write("variableStep chrom={} span=100\n".format(chrom))
            for summit, _, a_value, _ in tracks[chrom]:
                fout.write("{}\t{}\n".format(summit, a_value))

    file_name = os.path.join(ma.root_dir, "output_tracks", ma.output_prefix + "_P_values.wig")
    with open(file_name, "w") as fout:
        fout.write("track type=wiggle_0 name={}_-log10(P_value) visibility=full autoScale=on color=255,0,0 yLineMark=0 "
                   "yLineOnOff=on priority=10\n".format(ma.output_prefix))
        for chrom in tracks:
            if len(tracks[chrom]) == 0:
                continue
            fout.write("variableStep chrom={} span=100\n".format(chrom))
            for summit, _, _, p_value in tracks[chrom]:
                fout.write("{}\t{}\n".format(summit, -log10(p_value)))


def output_unbiased_peaks(ma, m_cutoff):
    peaks = defaultdict(list)
    for chrom in ma.peaks1.peaks:
        for peak in ma.peaks1.peaks[chrom]:
            if peak.type.endswith("unique"):
                peaks[chrom].append(peak)
    for chrom in ma.common_peaks.peaks:
        for peak in ma.common_peaks.peaks[chrom]:
            peaks[chrom].append(peak)
    for chrom in ma.peaks2.peaks:
        for peak in ma.peaks2.peaks[chrom]:
            if peak.type.endswith("unique"):
                peaks[chrom].append(peak)

    num = 0
    file_name = os.path.join(ma.root_dir, "output_filters", ma.output_prefix + "_unbiased_peaks.bed")
    with open(file_name, "w") as fout:
        for chrom in peaks:
            for peak in peaks[chrom]:
                if abs(peak.normed_m_value) < abs(m_cutoff):
                    num += 1
                    fout.write(
                        "{}\t{}\t{}\t{}\t{}\n".format(peak.chrom, peak.start, peak.end, peak.type, peak.normed_m_value))
    return num


def output_biased_peaks(ma, m_cutoff, p_cutoff):
    peaks = defaultdict(list)
    for chrom in ma.peaks1.peaks:
        for peak in ma.peaks1.peaks[chrom]:
            if peak.type.endswith("unique"):
                peaks[chrom].append(peak)
    for chrom in ma.common_peaks.peaks:
        for peak in ma.common_peaks.peaks[chrom]:
            peaks[chrom].append(peak)
    for chrom in ma.peaks2.peaks:
        for peak in ma.peaks2.peaks[chrom]:
            if peak.type.endswith("unique"):
                peaks[chrom].append(peak)

    num1, num2 = 0, 0
    file_name1 = os.path.join(ma.root_dir, "output_filters",
                              ma.peaks1.name + "_M_above_{}_biased_peaks.bed".format(abs(m_cutoff)))
    file_name2 = os.path.join(ma.root_dir, "output_filters",
                              ma.peaks2.name + "_M_below_-{}_biased_peaks.bed".format(abs(m_cutoff)))
    with open(file_name1, "w") as fout1, open(file_name2, "w") as fout2:
        for chrom in peaks:
            for peak in peaks[chrom]:
                if peak.p_value <= p_cutoff:
                    if peak.normed_m_value >= abs(m_cutoff):
                        num1 += 1
                        fout1.write(
                            "{}\t{}\t{}\t{}\t{}\n".format(peak.chrom, peak.start, peak.end, peak.type,
                                                          peak.normed_m_value))
                    elif peak.normed_m_value <= -abs(m_cutoff):
                        num2 += 1
                        fout2.write(
                            "{}\t{}\t{}\t{}\t{}\n".format(peak.chrom, peak.start, peak.end, peak.type,
                                                          peak.normed_m_value))

    return num1, num2
