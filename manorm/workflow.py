"""MAnorm workflow."""

import os
from manorm.logger import logger
from manorm.lib.io import *
from manorm.lib.peaks import *


def main(peaks_file1, peaks_file2, reads_file1, reads_file2, shift_size1, shift_size2, peak_width, distance_cutoff,
         random_times, overlap_dependent, m_cutoff, p_cutoff, output_all, output_name):
    if not os.path.isdir(output_name):
        os.mkdir(output_name)
    if not os.path.isdir(output_name + '/' + 'output_figures'):
        os.mkdir(output_name + '/' + 'output_figures')
    if not os.path.isdir(output_name + '/' + 'output_filters'):
        os.mkdir(output_name + '/' + 'output_filters')
    if not os.path.isdir(output_name + '/' + 'output_wig_files'):
        os.mkdir(output_name + '/' + 'output_wig_files')

    peaks_name1, peaks_name2 = os.path.basename(peaks_file1), os.path.basename(peaks_file2)
    reads_name1, reads_name2 = os.path.basename(reads_file1), os.path.basename(reads_file2)
    logger.info("# MAnorm Arguments:")
    logger.info("# Peaks file of sample1 = {}".format(peaks_name1))
    logger.info("# Peaks file of sample2 = {}".format(peaks_name2))
    logger.info("# Reads file of sample1 = {}".format(reads_file1))
    logger.info("# Reads file of sample2 = {}".format(reads_file2))
    logger.info("# Reads shift size of sample1 = {}".format(shift_size1))
    logger.info("# Reads shift size of sample2 = {}".format(shift_size2))
    logger.info("# Peak half width = {}".format(peak_width))
    logger.info("# Output folder = {}".format(output_name))

    peaks_name1 = os.path.splitext(peaks_name1)[0]
    peaks_name2 = os.path.splitext(peaks_name2)[0]
    reads_name1 = os.path.splitext(reads_name1)[0]
    reads_name2 = os.path.splitext(reads_name2)[0]

    logger.info("Step1: Loading input data...")
    peaks1, peaks2 = read_peaks(peaks_file1), read_peaks(peaks_file2)
    reads1, reads2 = read_reads(reads_file1, shift_size1), read_reads(reads_file2, shift_size2)

    logger.info("Step2: Classifying peaks by overlap...")
    peaks1_unique, peaks1_common, peaks2_unique, peaks2_common = get_common_peaks(peaks1, peaks2)
    logger.info(
        "{}: {}(unique) {}(common)".format(peaks_name1, get_peaks_size(peaks1_unique), get_peaks_size(peaks1_common)))
    logger.info(
        "{}: {}(unique) {}(common)".format(peaks_name2, get_peaks_size(peaks2_unique), get_peaks_size(peaks2_common)))
    logger.info("Performing overlap enrichment test...")
    fold_change = []
    for _ in range(random_times):
        tmp_peaks_random = randomize_peaks(peaks2)
        tmp_peaks_common = get_common_peaks(peaks1, tmp_peaks_random)[1]
        try:
            fold_change.append(1.0 * get_peaks_size(peaks1_common) / get_peaks_size(tmp_peaks_common))
        except ZeroDivisionError:
            fold_change.append(1.0 * get_peaks_size(peaks1_common))
    logger.info("Enrichment of peaks overlap: Fold Change: mean={0:f}, std={1:f}".format(np.array(fold_change).mean(),
                                                                                         np.array(fold_change).std()))

    logger.info("Step3: Merging common peaks...")
    merged_peaks, summit_dis = merge_common_peaks(peaks1_common, peaks2_common)
    logger.info("Merged common peaks: {}".format(get_peaks_size(merged_peaks)))

    logger.info("Step4: Calculating read density of each peak...")
    cal_peaks_read_density(peaks1, reads1, reads2, peak_width)
    cal_peaks_read_density(peaks2, reads1, reads2, peak_width)
    cal_peaks_read_density(merged_peaks, reads1, reads2, peak_width)

    logger.info("Step5: Fitting normalization model...")
    ma_fit = use_merged_peaks_fit_model(merged_peaks, summit_dis, distance_cutoff)
    if ma_fit[0] >= 0:
        logger.info("Normalization model: M = {0:f} * A + {1:f}".format(ma_fit[1], ma_fit[0]))
    else:
        logger.info("Normalization model: M = {0:f} * A - {1:f}".format(ma_fit[1], abs(ma_fit[0])))

    logger.info("Step6: Normalizing all peaks...")
    normalize_peaks(peaks1, ma_fit)
    normalize_peaks(peaks2, ma_fit)
    normalize_peaks(merged_peaks, ma_fit)

    logger.info("Step7: Output results...")
    if output_all:
        output_normalized_peaks(peaks1_unique, peaks1_common, output_name + '/' + peaks_name1 + '_MAvalues.xls',
                                reads_name1, reads_name2)
        output_normalized_peaks(peaks2_unique, peaks2_common, output_name + '/' + peaks_name2 + '_MAvalues.xls',
                                reads_name1, reads_name2)
    output_3set_normalized_peaks(peaks1_unique, merged_peaks, peaks2_unique,
                                 output_name + '/' + output_name + '_all_peaks_MAvalues.xls',
                                 peaks_name1, peaks_name2, reads_name1, reads_name2)

    draw_figs_to_show_data(output_name + '/' + 'output_figures', peaks1_unique, peaks2_unique, merged_peaks,
                           peaks_name1, peaks_name2, ma_fit, reads_name1, reads_name2)
    output_peaks_mvalue_2wig_file(output_name + '/' + 'output_wig_files', peaks1_unique, peaks2_unique, merged_peaks,
                                  output_name)
    unbiased_mvalue = m_cutoff
    output_unbiased_peaks(output_name + '/' + 'output_filters', peaks1_unique, peaks2_unique, merged_peaks,
                          unbiased_mvalue, overlap_dependent)
    output_biased_peaks(output_name + '/' + 'output_filters', peaks1_unique, peaks2_unique, merged_peaks, m_cutoff,
                        p_cutoff, overlap_dependent)


if __name__ == '__main__':
    pass
