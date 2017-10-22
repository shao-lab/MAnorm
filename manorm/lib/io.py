"""IO for MAnorm."""

import numpy as np
from numpy.ma import log2, log10
import matplotlib

matplotlib.use('Agg')
from matplotlib import pyplot as plt
from manorm.logger import logger
from manorm.lib.peaks import Peak, get_peaks_mavalues, get_peaks_normed_mavalues, get_peaks_pvalues, _add_peaks, \
    _sort_peaks_list


def _get_reads_position(reads_fp, shift):
    position = {}
    with open(reads_fp) as fin:
        for line in fin:
            fields = line.split('\t')
            chrm, start, end, strand = fields[0].strip(), int(fields[1]), int(fields[2]), fields[5].strip()
            pos = start + shift if strand is '+' else end - shift
            try:
                position[chrm].append(pos)
            except KeyError:
                position[chrm] = []
                position[chrm].append(pos)
    return {chrom: sorted(position[chrom]) for chrom in position.keys()}


def _get_read_length(reads_file):
    with open(reads_file) as fin:
        for line in fin:
            fields = line.split('\t')
            return int(fields[2]) - int(fields[1])


def read_reads(reads_file, shift):
    return _get_reads_position(reads_file, shift)


def _read_peaks(peaks_file):
    peaks = {}
    with open(peaks_file) as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            chrom = fields[0].strip()
            try:
                peak = Peak(chrom, int(fields[1]), int(fields[2]), int(fields[3].strip()))
            except:
                peak = Peak(chrom, int(fields[1]), int(fields[2]))
            try:
                peaks[chrom].append(peak)
            except KeyError:
                peaks[chrom] = []
                peaks[chrom].append(peak)
    return peaks


def _read_macs_xls_peaks(peaks_file):
    peaks = {}
    with open(peaks_file) as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            chrom = fields[0].strip()
            try:
                peak = Peak(chrom, int(fields[1]), int(fields[2]), int(fields[4]))
            except:
                continue
            try:
                peaks[chrom].append(peak)
            except KeyError:
                peaks[chrom] = []
                peaks[chrom].append(peak)
    return peaks


def read_peaks(peaks_file):
    try:
        return _read_peaks(peaks_file)
    except:
        return _read_macs_xls_peaks(peaks_file)


def output_normalized_peaks(peaks_unique, peaks_common, file_name, reads1_name, reads2_name):
    fout = open(file_name, 'w')
    header = '\t'.join(['chr', 'start', 'end', 'summit', 'M_value', 'A_value', 'P_value', 'Peak_Group',
                        'normalized_read_density_in_%s' % reads1_name, 'normalized_read_density_in_%s\n' % reads2_name])
    fout.write(header)
    for chrom in peaks_unique.keys():
        for peak in peaks_unique[chrom]:
            cnt = (peak.chrm, peak.start, peak.end, peak.summit - peak.start,
                   peak.normed_mvalue, peak.normed_avalue, str(peak.pvalue), 'unique',
                   peak.normed_read_density1, peak.normed_read_density2,)
            fout.write('\t'.join(['%s', '%d', '%d', '%d', '%f', '%f', '%s', '%s', '%f', '%f']) % cnt + '\n')
    for chrom in peaks_common.keys():
        for peak in peaks_common[chrom]:
            cnt = (peak.chrm, peak.start, peak.end, peak.summit - peak.start,
                   peak.normed_mvalue, peak.normed_avalue, str(peak.pvalue), 'common',
                   peak.normed_read_density1, peak.normed_read_density2)
            fout.write('\t'.join(['%s', '%d', '%d', '%d', '%f', '%f', '%s', '%s', '%f', '%f']) % cnt + '\n')
    fout.close()


def output_3set_normalized_peaks(peaks1_unique, merged_peaks, peaks2_unique, file_name, peaks1_name, peaks2_name,
                                 reads1_name, reads2_name):
    fout = open(file_name, 'w')
    header = '\t'.join(['chr', 'start', 'end', 'summit', 'M_value', 'A_value', 'P_value', 'Peak_Group',
                        'normalized_read_density_in_%s' % reads1_name, 'normalized_read_density_in_%s\n' % reads2_name])
    fout.write(header)
    for chrom in peaks1_unique.keys():
        for peak in peaks1_unique[chrom]:
            cnt = (peak.chrm, peak.start, peak.end, peak.summit - peak.start,
                   peak.normed_mvalue, peak.normed_avalue, str(peak.pvalue), '%s_unique' % peaks1_name,
                   peak.normed_read_density1, peak.normed_read_density2)
            fout.write('\t'.join(['%s', '%d', '%d', '%d', '%f', '%f', '%s', '%s', '%f', '%f']) % cnt + '\n')
    for chrom in merged_peaks.keys():
        for peak in merged_peaks[chrom]:
            cnt = (peak.chrm, peak.start, peak.end, peak.summit - peak.start,
                   peak.normed_mvalue, peak.normed_avalue, str(peak.pvalue), 'merged_common_peak',
                   peak.normed_read_density1, peak.normed_read_density2)
            fout.write('\t'.join(['%s', '%d', '%d', '%d', '%f', '%f', '%s', '%s', '%f', '%f']) % cnt + '\n')
    for chrom in peaks2_unique.keys():
        for peak in peaks2_unique[chrom]:
            cnt = (peak.chrm, peak.start, peak.end, peak.summit - peak.start,
                   peak.normed_mvalue, peak.normed_avalue, str(peak.pvalue), '%s_unique' % peaks2_name,
                   peak.normed_read_density1, peak.normed_read_density2)
            fout.write('\t'.join(['%s', '%d', '%d', '%d', '%f', '%f', '%s', '%s', '%f', '%f']) % cnt + '\n')
    fout.close()


def draw_figs_to_show_data(output_dir, peaks1_unique, peaks2_unique, merged_peaks, peaks1_name, peaks2_name, ma_fit,
                           reads1_name, reads2_name):
    peaks_3set = [peaks1_unique, peaks2_unique, merged_peaks]
    peaks1_name = ' '.join([peaks1_name, 'unique'])
    peaks2_name = ' '.join([peaks2_name, 'unique'])
    merged_pks_name = 'merged common peaks'
    peaks_names = [peaks1_name, peaks2_name, merged_pks_name]
    colors = 'bgr'
    a_max = 0
    a_min = 10000
    plt.figure(1).set_size_inches(16, 12)
    for (idx, peaks) in enumerate(peaks_3set):
        mvalues, avalues = get_peaks_mavalues(peaks)
        if len(avalues) != 0:
            a_max = max(max(avalues), a_max)
            a_min = min(min(avalues), a_min)
        plt.scatter(avalues, mvalues, s=10, c=colors[idx])
    plt.xlabel('A value')
    plt.ylabel('M value')
    plt.grid(axis='y')
    plt.legend(peaks_names, loc='best')
    plt.title('before rescale')
    x = np.arange(a_min, a_max, 0.01)
    y = ma_fit[1] * x + ma_fit[0]
    plt.plot(x, y, '-', color='k')
    plt.savefig(output_dir + '/' + 'before_rescale.png')

    plt.figure(2).set_size_inches(16, 12)
    rd_min = 1000000
    rd_max = 0
    reads_density1, reads_density2 = [], []
    for chrom in merged_peaks.keys():
        for peak in merged_peaks[chrom]:
            reads_density1.append(peak.read_density1), reads_density2.append(peak.read_density2)
    rd_max = max(max(log2(reads_density1)), rd_max)
    rd_min = min(min(log2(reads_density1)), rd_min)
    plt.scatter(log2(reads_density1), log2(reads_density2), s=10, c='r', label=merged_pks_name, alpha=0.5)
    plt.xlabel(' log2 read density' + ' by ' + '"' + reads1_name + '" reads')
    plt.ylabel(' log2 read density' + ' by ' + '"' + reads2_name + '" reads')
    plt.grid(axis='y')
    plt.legend(loc='upper left')
    plt.title('Fitting Model via common peaks')
    rx = np.arange(rd_min, rd_max, 0.01)
    ry = (2 - ma_fit[1]) * rx / (2 + ma_fit[1]) - 2 * ma_fit[0] / (2 + ma_fit[1])
    plt.plot(rx, ry, '-', color='k')
    plt.savefig(output_dir + '/' + 'log2_read_density.png')

    # plot the MA plot after rescale
    plt.figure(3).set_size_inches(16, 12)
    for (idx, peaks) in enumerate(peaks_3set):
        normed_mvalues, normed_avalues = get_peaks_normed_mavalues(peaks)
        plt.scatter(normed_avalues, normed_mvalues, s=10, c=colors[idx])
    plt.xlabel('A value')
    plt.ylabel('M value')
    plt.grid(axis='y')
    plt.legend(peaks_names, loc='best')
    plt.title('after rescale')
    plt.savefig(output_dir + '/' + 'after_rescale.png')

    plt.figure(4).set_size_inches(16, 12)
    for (idx, peaks) in enumerate(peaks_3set):
        normed_mvalues, normed_avalues = get_peaks_normed_mavalues(peaks)
        colors = -log10(get_peaks_pvalues(peaks))
        for i, c in enumerate(colors):
            if c > 50:
                colors[i] = 50
        plt.scatter(normed_avalues, normed_mvalues, s=10, c=colors, cmap='jet')
    plt.colorbar()
    plt.grid(axis='y')
    plt.xlabel('A value')
    plt.ylabel('M value')
    plt.title('-log10(P-value)')
    plt.savefig(output_dir + '/' + '-log10_P-value.png')
    plt.close()


def output_peaks_mvalue_2wig_file(output_dir, peaks1_unique, peaks2_unique, merged_peaks, comparison_name):
    peaks = _add_peaks(_add_peaks(peaks1_unique, merged_peaks), peaks2_unique)
    fout = open(output_dir +'/'+ '_'.join([comparison_name, 'peaks_Mvalues.wig']), 'w')
    fout.write('browser position chr11:5220000-5330000\n')
    fout.write('track type=wiggle_0 name=%s' % comparison_name + ' visibility=full autoScale=on color=255,0,0 ' +
               ' yLineMark=0 yLineOnOff=on priority=10\n')
    for chrom in peaks.keys():
        fout.write('variableStep chrom=' + chrom + ' span=100\n')
        peaks_chr = peaks[chrom]
        sorted_pks_chr = _sort_peaks_list(peaks_chr, 'summit')
        # write sorted peak summit and m-value to file
        [fout.write('\t'.join(['%d' % peak.summit, '%s\n' % str(peak.normed_mvalue)])) for peak in sorted_pks_chr]
    fout.close()

    fout = open(output_dir +'/'+ '_'.join([comparison_name, 'peaks_Pvalues.wig']), 'w')
    fout.write('browser position chr11:5220000-5330000\n')
    fout.write('track type=wiggle_0 name=%s(-log10(p-value))' % comparison_name +
               ' visibility=full autoScale=on color=255,0,0 ' + ' yLineMark=0 yLineOnOff=on priority=10\n')
    for chrom in peaks.keys():
        fout.write('variableStep chrom=' + chrom + ' span=100\n')
        peaks_chr = peaks[chrom]
        sorted_pks_chr = _sort_peaks_list(peaks_chr, 'summit')
        # write sorted peak summit and m-value to file
        [fout.write('\t'.join(['%d' % peak.summit, '%s\n' % str(-log10(peak.pvalue))])) for peak in sorted_pks_chr]
    fout.close()


def output_unbiased_peaks(output_dir, peaks1_unique, peaks2_unique, merged_peaks, unbiased_mvalue, overlap_dependent):
    if not overlap_dependent:
        pks = _add_peaks(_add_peaks(peaks1_unique, merged_peaks), peaks2_unique)
        name = 'all_peaks'
    else:
        pks = merged_peaks
        name = 'merged_common_peaks'

    file_bed = open(output_dir + '/' + 'unbiased_peaks_of_%s' % name + '.bed', 'w')
    # file_bed.write(bed_peak_header)
    i = 0
    for key in pks.keys():
        for pk in pks[key]:
            if abs(pk.normed_mvalue) < unbiased_mvalue:
                i += 1
                line = '\t'.join([pk.chrm, '%d' % pk.start, '%d' % pk.end, 'from_%s_%d' % (name, i),
                                  '%s\n' % str(pk.normed_mvalue)])
                file_bed.write(line)
    logger.info('Unbiased peaks: {}'.format(i))
    file_bed.close()


def output_biased_peaks(output_dir, peaks1_unique, peaks2_unique, merged_pks, biased_mvalue, biased_pvalue,
                        overlap_dependent):
    if not overlap_dependent:
        peaks = _add_peaks(_add_peaks(peaks1_unique, merged_pks), peaks2_unique)
        name = 'all_peaks'
    else:
        peaks = _add_peaks(peaks1_unique, peaks2_unique)
        name = 'unique_peaks'

    file_bed_over = open(output_dir + '/' + 'M_over_%.2f_biased_peaks_of_%s' % (biased_mvalue, name) + '.bed', 'w')
    # file_bed_over.write(bed_peak_header)
    file_bed_less = open(output_dir + '/' + 'M_less_-%.2f_biased_peaks_of_%s' % (biased_mvalue, name) + '.bed', 'w')
    # file_bed_less.write(bed_peak_header)
    i, j = 0, 0
    for chrom in peaks.keys():
        for peak in peaks[chrom]:
            if peak.pvalue < biased_pvalue:
                if peak.normed_mvalue > biased_mvalue:
                    i += 1
                    line = '\t'.join([peak.chrm, '%d' % peak.start, '%d' % peak.end, 'from_%s_%d' % (name, i),
                                      '%s\n' % str(peak.normed_mvalue)])
                    file_bed_over.write(line)
                if peak.normed_mvalue < -biased_mvalue:
                    j += 1
                    line = '\t'.join([peak.chrm, '%d' % peak.start, '%d' % peak.end, 'from_%s_%d' % (name, j),
                                      '%s\n' % str(peak.normed_mvalue)])
                    file_bed_less.write(line)
    logger.info('Biased peaks: {}'.format(i + j))
    file_bed_over.close()
    file_bed_less.close()
