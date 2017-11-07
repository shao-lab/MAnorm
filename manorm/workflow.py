"""MAnorm workflow."""

import os
from manorm.core import MAnorm
from manorm.logger import logger
from manorm.io import *
from manorm.plot import ma_plot


def main(peaks_file1, peaks_file2, reads_file1, reads_file2, shift_size1, shift_size2, peak_width, summit_dis_cutoff,
         n_random, m_cutoff, p_cutoff, full_output, name1, name2, output_name):
    root_dir = output_name
    output_prefix = os.path.basename(output_name)
    logger.info("# MAnorm Arguments:")
    if name1 and name2:
        logger.info("# Name of sample1 = {}".format(name1))
        logger.info("# Name of sample2 = {}".format(name2))
    logger.info("# Peaks file of sample1 = {}".format(os.path.basename(peaks_file1)))
    logger.info("# Peaks file of sample2 = {}".format(os.path.basename(peaks_file2)))
    logger.info("# Reads file of sample1 = {}".format(os.path.basename(reads_file1)))
    logger.info("# Reads file of sample2 = {}".format(os.path.basename(reads_file2)))
    logger.info("# Reads shift size of sample1 = {}".format(shift_size1))
    logger.info("# Reads shift size of sample2 = {}".format(shift_size2))
    logger.info("# Peak width to extend from summit = {}".format(peak_width))
    logger.info("# Summit-to-summit distance cutoff = {}".format(summit_dis_cutoff))
    logger.info("# Output directory = {}".format(root_dir))

    mk_dir(root_dir)
    ma = MAnorm(name1, name2, root_dir, output_prefix)

    logger.info("Step1: Loading input data")
    ma.load_peaks(peaks_file1, peaks_file2)
    ma.load_reads(reads_file1, reads_file2, shift_size1, shift_size2)

    logger.info("Step2: Classifying peaks by overlap")
    ma.classify_peaks_by_overlap()
    ma.merge_common_peaks()
    unique_num1 = ma.peaks1.count_peak_type("unique")
    common_num1 = ma.peaks1.count_peak_type("common")
    unique_num2 = ma.peaks2.count_peak_type("unique")
    common_num2 = ma.peaks2.count_peak_type("common")
    logger.info("{}: {}(unique) {}(common)".format(ma.peaks1.name, unique_num1, common_num1))
    logger.info("{}: {}(unique) {}(common)".format(ma.peaks2.name, unique_num2, common_num2))
    logger.info("Merged common peaks: {}".format(ma.common_peaks.size))

    logger.info("Step3: Performing overlap enrichment test")
    random_overlap_num, fold_change = ma.overlap_enrichment_test(n_random)
    logger.info("Overlapping peaks in random simulation: mean:{} std:{}".format(random_overlap_num.mean(),
                                                                                random_overlap_num.std()))
    logger.info("Fold Change: {}".format(fold_change.mean()))

    logger.info("Step4: Calculating the read densities")
    ma.cal_read_density(peak_width)

    logger.info("Step5: Fitting normalization model")
    ma.fit_model(summit_dis_cutoff)
    ma_params = ma.ma_model
    if ma_params[0] >= 0:
        logger.info("M-A model: M = {0:f} * A + {1:f}".format(ma_params[1], ma_params[0]))
    else:
        logger.info("M-A model: M = {0:f} * A - {1:f}".format(ma_params[1], abs(ma_params[0])))

    logger.info("Step6: Normalizing all peaks")
    ma.normalize()

    logger.info("Step7: Output results")
    if full_output:
        output_original_peaks(ma)
    output_all_peaks(ma)
    output_wiggle_track(ma)
    unbiased_num = output_unbiased_peaks(ma, m_cutoff)
    logger.info("{} peaks with |M_value|<{} are filtered as unbiased peaks".format(unbiased_num, abs(m_cutoff)))
    biased_num1, biased_num2 = output_biased_peaks(ma, m_cutoff, p_cutoff)
    logger.info("{} peaks with M_value>={} are filtered as sample1-biased peaks".format(biased_num1, abs(m_cutoff)))
    logger.info("{} peaks with M_value<=-{} are filtered as sample2-biased peaks".format(biased_num2, abs(m_cutoff)))
    ma_plot(ma)
    logger.info("Finished!")
