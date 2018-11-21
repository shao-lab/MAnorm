"""This module implements the MAnorm application."""
from __future__ import absolute_import

import logging
import sys
import time
from os import path

import manorm
from manorm import peak, read
from manorm.logger import setup_logger
from manorm.model import MAmodel

logger = logging.getLogger(__name__)


class MAnorm(object):
    """The abstract application class for MAnorm."""

    def __init__(self, args):
        """Initialize MAnorm application.

        :param args: The ``argparse.Namespace`` object which contains arguments for running MAnorm.
        """
        self.start_time = time.time()
        self.end_time = None
        self.version = manorm.__version__
        self.args = args
        # set default value for summit-dis cutoff
        if self.args.summit_dis is None:
            self.args.summit_dis = self.args.window_size // 4
        self.name1 = args.name1 if args.name1 else path.splitext(path.basename(args.peaks_file1))[0]
        self.name2 = args.name2 if args.name2 else path.splitext(path.basename(args.peaks_file2))[0]
        self.output_prefix = self.name1 + "_vs_" + self.name2
        self.output_dir = args.output_dir
        self.ma_model = None

    def print_args(self):
        logger.info("")
        logger.info("==== Parameters =====")
        logger.info("Sample 1 name: {}".format(self.name1))
        logger.info("Sample 1 peaks file: {}".format(self.args.peaks_file1))
        logger.info("Sample 1 reads file: {}".format(self.args.reads_file1))
        if not self.args.paired:
            logger.info("Sample 1 reads shift size: {}".format(self.args.shift_size1))
        logger.info("Sample 2 name: {}".format(self.name2))
        logger.info("Sample 2 peaks file: {}".format(self.args.peaks_file2))
        logger.info("Sample 2 reads file: {}".format(self.args.reads_file2))
        if not self.args.paired:
            logger.info("Sample 2 reads shift size: {}".format(self.args.shift_size2))
            logger.info("Reads mode: single-end")
        else:
            logger.info("Reads mode: paired-end")
        logger.info("Window size to count reads: {}".format(self.args.window_size))
        logger.info("Summit distance cutoff: {}".format(self.args.summit_dis))
        logger.info("Number of random simulation: {}".format(self.args.n_random))
        logger.info("M-value cutoff: {}".format(self.args.m_cutoff))
        logger.info("P-value cutoff: {}".format(self.args.p_cutoff))
        logger.info("Output directory: {}".format(self.args.output_dir))

    def load_input_data(self):
        """Load input data and initiate the MA model."""
        logger.debug("Loading peaks of sample 1")
        peaks1 = peak.load_peaks(path=self.args.peaks_file1, format=self.args.peaks_format, name=self.name1)
        logger.debug("Loaded peaks: {}".format(peaks1.size))

        logger.debug("Loading peaks of sample 2")
        peaks2 = peak.load_peaks(path=self.args.peaks_file2, format=self.args.peaks_format, name=self.name2)
        logger.debug("Loaded peaks: {}".format(peaks2.size))

        logger.debug("Loading reads of sample 1")
        reads1 = read.load_reads(path=self.args.reads_file1, paired=self.args.paired, shift=self.args.shift_size1)
        logger.debug("Loaded reads: {:,}".format(reads1.size))

        logger.debug("Loading reads of sample 2")
        reads2 = read.load_reads(path=self.args.reads_file2, paired=self.args.paired, shift=self.args.shift_size2)
        logger.debug("Loaded reads: {:,}".format(reads2.size))

        logging.debug("Initiating MA model")
        self.ma_model = MAmodel(peaks1=peaks1, peaks2=peaks2, reads1=reads1, reads2=reads2)
        logging.debug("MA model initiated")

    def run(self):
        setup_logger(self.args.verbose)
        logger.info("Running MAnorm v{}".format(self.version))
        logger.info("Command-line arguments: {}".format(' '.join(sys.argv[1:])))
        self.print_args()
        logger.info("")
        logger.info("==== Running ====")

        logger.info("Loading input data")
        self.load_input_data()

        logger.info("Classifying peaks by overlap")
        self.ma_model.classify_peaks_by_overlap()
        self.ma_model.merge_common_peaks()

        if self.args.n_random > 0:
            logger.info("Testing the enrichment of peak overlap")
            # TODO: overlap enrichment test
            self.ma_model.overlap_enrichment_test(self.args.n_random)

        # TODO: manorm-pipeline

        logger.info("Calculating the read densities")
        self.ma_model.cal_read_density(self.args.peak_width)

        logger.info("Fitting normalization model")
        self.ma_model.fit_model(self.args.summit_dis_cutoff)
        logger.info("Normalizing all peaks")
        self.ma_model.normalize()

        logger.info("Output results")
        self.output()


    def output(self):
        # output_original_peaks(ma)
        # output_all_peaks(ma)
        # output_wiggle_track(ma)
        # output_unbiased_peaks(ma, m_cutoff)
        # output_biased_peaks(ma, m_cutoff, p_cutoff)
        # ma_plot(ma)
        pass

    def report(self):
        logger.info("==== Statistics ====")
        # logger.info("Total reads of sample 1: {:,}".format(len(self.reads1)))
        # logger.info("Total reads of sample 2: {:,}".format(len(self.reads2)))
        # logger.info("Total peaks of sample 1: {} (unique: {} common: {})".format(len(self.peaks1), 100, 200))
        # logger.info("Total peaks of sample 2: {} (unique: {} common: {})".format(len(self.peaks2), 200, 100))
        # logger.info("Merged common peaks: {}".format(len(self.peaks_common)))
        # logger.info("Overlapping peaks in random simulation: mean:{} std:{}".format(random_overlap_num.mean(),
        #                                                                             random_overlap_num.std()))
        # logger.info("Fold Change: {}".format(fold_change.mean()))
        # if ma_params[0] >= 0:
        #     logger.info("M-A model: M = {0:f} * A + {1:f}".format(ma_params[1], ma_params[0]))
        # else:
        #     logger.info("M-A model: M = {0:f} * A - {1:f}".format(ma_params[1], abs(ma_params[0])))
        # logger.info("{} peaks with |M_value|<{} are filtered as unbiased peaks".format(unbiased_num, abs(m_cutoff)))
        # logger.info("{} peaks with M_value>={} are filtered as sample1-biased peaks".format(biased_num1, abs(m_cutoff)))
        # logger.info(
        #     "{} peaks with M_value<=-{} are filtered as sample2-biased peaks".format(biased_num2, abs(m_cutoff)))
        pass

    def exit(self):
        pass
