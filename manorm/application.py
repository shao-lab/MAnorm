# -*- coding: utf-8 -*-

"""
manorm.application
~~~~~~~~~~~~~~~~~~

MAnorm application/pipeline.
"""

from __future__ import absolute_import, division

import logging
import sys
import time
from os import path

import numpy as np

import manorm
from manorm import peak, read
from manorm.compat import range
from manorm.logging import setup_logger
from manorm.model import MAmodel
from manorm.peak.utils import generate_random_peaks, overlap_on_single_chr

logger = logging.getLogger(__name__)


class MAnorm(object):
    """The application class for MAnorm."""

    def __init__(self, peak_file1, peak_file2, read_file1, read_file2, peak_format='bed', read_format='bed',
                 name1=None, name2=None, shift_size1=100, shift_size2=100, paired=False, window_size=2000,
                 summit_dis=None, n_random=10, m_cutoff=1.0, p_cutoff=0.01, output_all=False, output_dir=None,
                 verbose=False, *args, **kwargs):
        """Initialize MAnorm application.

        :param peak_file1: Path of sample 1 peak file.
        :param peak_file2: Path of sample 2 peak file.
        :param read_file1: Path of sample 1 read file.
        :param read_file2: Path of sample 2 read file.
        :param peak_format: Format of peak files.
        :param read_format: Format of read files.
        :param name1: Name of sample 1.
        :param name2: Name of sample 2.
        :param shift_size1: Read shift size of sample 1, only effective in single-end mode.
        :param shift_size2: Read shift size of sample 2, only effective in single-end mode.
        :param paired: Paired-end sequencing or not.
        :param window_size: Window size when counting reads.
        :param summit_dis: Summit distance cutoff for overlapping common peaks when fitting normalization model.
        :param n_random: Number of random simulations to test the enrichment of peak overlapping.
        :param m_cutoff: M value cutoff to define biased (differential binding) peaks.
        :param p_cutoff: P value cutoff to define biased (differential binding) peaks.
        :param output_all: Write two extra output files containing the results of original (unmerged) peaks.
        :param output_dir: Output directory.
        :param verbose: Enable verbose log messages or not.
        :param args: Other positional arguments.
        :param kwargs: Other keyword arguments.
        """
        self.start_time = time.time()
        self.end_time = None
        self.version = manorm.__version__
        self.peak_file1 = path.abspath(peak_file1)
        self.peak_file2 = path.abspath(peak_file2)
        self.peak_format = peak_format
        self.read_file1 = path.abspath(read_file1)
        self.read_file2 = path.abspath(read_file2)
        self.read_format = read_format
        self.shift_size1 = shift_size1
        self.shift_size2 = shift_size2
        self.paired = paired
        self.window_size = window_size
        self.summit_dis = window_size // 4 if summit_dis is None else summit_dis
        self.n_random = n_random
        self.n_overlap_rand = None
        self.m_cutoff = m_cutoff
        self.p_cutoff = p_cutoff
        self.name1 = name1 if name1 else path.splitext(path.basename(peak_file1))[0]
        self.name2 = name2 if name2 else path.splitext(path.basename(peak_file2))[0]
        self.output_prefix = self.name1 + "_vs_" + self.name2
        self.output_all = output_all
        self.output_dir = path.abspath('') if output_dir is None else path.abspath(output_dir)
        self.verbose = verbose
        self.ma_model = None
        self.args = args
        self.kwargs = kwargs

    def print_args(self):
        """Log parameters that MAnorm uses."""
        logger.info("==== INFO ====")
        logger.info("Sample 1 name = {}".format(self.name1))
        logger.info("Sample 1 peak file = {} [{}]".format(self.peak_file1, self.peak_format))
        logger.info("Sample 1 read file = {} [{}]".format(self.read_file1, self.read_format))
        if not self.paired:
            logger.info("Sample 1 read shift size = {}".format(self.shift_size1))
        logger.info("Sample 2 name = {}".format(self.name2))
        logger.info("Sample 2 peak file = {} [{}]".format(self.peak_file2, self.peak_format))
        logger.info("Sample 2 read file = {} [{}]".format(self.read_file2, self.read_format))
        if not self.paired:
            logger.info("Sample 2 read shift size = {}".format(self.shift_size2))
        if self.paired:
            logger.info("Paired-end mode: on")
        else:
            logger.info("Paired-end mode: off")
        logger.info("Window size = {}".format(self.window_size))
        logger.info("Summit distance cutoff = {}".format(self.summit_dis))
        logger.info("Number of random simulation = {}".format(self.n_random))
        logger.info("M-value cutoff = {}".format(self.m_cutoff))
        logger.info("P-value cutoff = {}".format(self.p_cutoff))
        logger.info("Output directory = {}".format(self.output_dir))

    def load_input_data(self):
        """Read input data."""
        logger.info("Loading peaks of sample 1")
        peaks1 = peak.load_peaks(path=self.peak_file1, format=self.peak_format, name=self.name1)
        logger.info("Loading peaks of sample 2")
        peaks2 = peak.load_peaks(path=self.peak_file2, format=self.peak_format, name=self.name2)
        logger.info("Loading reads of sample 1")
        reads1 = read.load_reads(path=self.read_file1, format=self.read_format, paired=self.paired,
                                 shift=self.shift_size1, name=self.name1)
        logger.info("Loading reads of sample 2")
        reads2 = read.load_reads(path=self.read_file2, format=self.read_format, paired=self.paired,
                                 shift=self.shift_size2, name=self.name2)
        return peaks1, peaks2, reads1, reads2

    def test_overlap_enrich(self, peaks1, peaks2):
        """Given two sets of peaks, test the enrichment of overlap compared to random control peak sets."""
        n_overlap_rand = []
        for _ in range(self.n_random):
            peak_rand = generate_random_peaks(peaks2)
            tmp_count = 0
            for chrom in set(peaks1.chroms) | set(peak_rand.chroms):
                flag_overlap, _ = overlap_on_single_chr(peaks1.fetch(chrom), peak_rand.fetch(chrom))
                tmp_count += flag_overlap.sum()
            if tmp_count == 0:
                tmp_count = 1
            n_overlap_rand.append(tmp_count)
        n_overlap_rand = np.array(n_overlap_rand)
        return n_overlap_rand

    def run(self):
        """Run MAnorm pipeline."""
        setup_logger(self.verbose)
        logger.info("Running MAnorm v{}".format(self.version))
        logger.info("Command-line arguments: {}".format(' '.join(sys.argv[1:])))
        self.print_args()

        logger.info("==== Running ====")
        logger.info("Step 1: Loading input data")
        peaks1, peaks2, reads1, reads2 = self.load_input_data()

        logger.info("Step 2: Processing peaks")
        self.ma_model = MAmodel(peaks1, peaks2, reads1, reads2)
        self.ma_model.process_peaks()

        logger.info("Step 3: Testing the enrichment of peak overlap")
        if self.n_random > 0:
            self.n_overlap_rand = self.test_overlap_enrich(peaks1, peaks2)
        else:
            logger.info("Skipped")

        logger.info("Step 4: Fitting M-A normalization model on common peaks")
        self.ma_model.fit_model(window_size=self.window_size, summit_dis=self.summit_dis)

        logger.info("Step 5: Normalizing all peaks")
        self.ma_model.normalize()

        logger.info("Step 6: Write output files")
        self.output()

        self.report()

    def output(self):
        """Write output files."""
        # output_original_peaks(ma)
        # output_all_peaks(ma)
        # output_wiggle_track(ma)
        # output_unbiased_peaks(ma, m_cutoff)
        # output_biased_peaks(ma, m_cutoff, p_cutoff)
        # ma_plot(ma)
        pass

    def report(self):
        """Report statistics."""
        logger.info("==== Stats ====")
        if self.paired:
            logger.info("Total read pairs of sample 1: {:,}".format(self.ma_model.reads1.size))
            logger.info("Total read pairs of sample 2: {:,}".format(self.ma_model.reads2.size))
        else:
            logger.info("Total reads of sample 1: {:,}".format(self.ma_model.reads1.size))
            logger.info("Total reads of sample 2: {:,}".format(self.ma_model.reads2.size))
        logger.info("Total peaks of sample 1: {} (unique: {} common: {})".format(
            self.ma_model.peaks1.size, self.ma_model.peaks1.n_unique, self.ma_model.peaks1.n_common))
        logger.info("Total peaks of sample 2: {} (unique: {} common: {})".format(
            self.ma_model.peaks2.size, self.ma_model.peaks2.n_unique, self.ma_model.peaks2.n_common))
        logger.info("Number of merged common peaks: {}".format(self.ma_model.peaks_merged.size))

        logger.info("Number of overlapping peaks in random: mean:{:.2f} std:{:.2f}".format(
            self.n_overlap_rand.mean(), self.n_overlap_rand.std()))
        logger.info("Fold change of overlapping peaks compared to random: {:.2f}".format(
            self.ma_model.peaks1.n_common / self.n_overlap_rand.mean()))
        if self.ma_model.ma_params[0] >= 0:
            logger.info("M-A model: M = {:f} * A + {:f}".format(self.ma_model.ma_params[1], self.ma_model.ma_params[0]))
        else:
            logger.info("M-A model: M = {:f} * A - {:f}".format(self.ma_model.ma_params[1], self.ma_model.ma_params[0]))
        # logger.info("{} peaks with |M_value|<{} are filtered as unbiased peaks".format(unbiased_num, abs(m_cutoff)))
        # logger.info("{} peaks with M_value>={} are filtered as sample1-biased peaks".format(biased_num1, abs(m_cutoff)))
        # logger.info(
        #     "{} peaks with M_value<=-{} are filtered as sample2-biased peaks".format(biased_num2, abs(m_cutoff)))
