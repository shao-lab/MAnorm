"""
manorm.cli
----------

Command line interface of MAnorm.
"""

import argparse
import logging
import os
from textwrap import dedent

from manorm import __version__
from manorm.io import mk_dir, write_all_peaks, write_original_peaks, \
    write_biased_peaks, write_wiggle_track
from manorm.logging import setup_logger
from manorm.model import MAmodel
from manorm.plot import plt_figures
from manorm.read import READ_FORMATS, load_reads
from manorm.region import REGION_FORMATS, load_manorm_peaks
from manorm.region.utils import random_peak_overlap, count_common_peaks, \
    count_unique_peaks

logger = logging.getLogger(__name__)


def _existed_file(path):
    """Check whether a passed argument is an existed file."""
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"file not found: {path}")
    return path


def _pos_int(value):
    """Check whether a passed argument is a positive integer."""
    try:
        value_int = int(value)
        if value_int <= 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"invalid positive int value: {value!r}")
    return value_int


def configure_parser():
    """Configure the arguments parser for MAnorm."""
    description = dedent("""
    MAnorm -- A robust model for quantitative comparison of ChIP-seq data sets
   
    Citation:
      Shao Z, Zhang Y, Yuan GC, Orkin SH, Waxman DJ. MAnorm: a robust model 
      for quantitative comparison of ChIP-Seq data sets. Genome biology. 
      2012 Mar 16;13(3):R16. https://doi.org/10.1186/gb-2012-13-3-r16
    """)

    epilog = dedent("""    
    See also:
      Documentation: https://manorm.readthedocs.io
      Source code: https://github.com/shao-lab/MAnorm
      Bug reports: https://github.com/shao-lab/MAnorm/issues
    """)

    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "-v", "--version", action="version", version=f"MAnorm {__version__}")

    parser_input = parser.add_argument_group("Input Options")
    parser_input.add_argument(
        "--p1", "--peak1", metavar="FILE", dest="peak_file1", required=True,
        type=_existed_file, help="Peak file of sample 1.")
    parser_input.add_argument(
        "--p2", "--peak2", metavar="FILE", dest="peak_file2", required=True,
        type=_existed_file, help="Peak file of sample 2.")
    parser_input.add_argument(
        "--pf", "--peak-format", metavar="FORMAT", dest="peak_format",
        choices=REGION_FORMATS, default="bed",
        help=f"Format of the peak files. Support {REGION_FORMATS}. "
             f"Default: bed")
    parser_input.add_argument(
        "--r1", "--read1", metavar="FILE", dest="read_file1", required=True,
        type=_existed_file, help="Read file of sample 1.")
    parser_input.add_argument(
        "--r2", "--read2", metavar="FILE", dest="read_file2", required=True,
        type=_existed_file, help="Read file of sample 2.")
    parser_input.add_argument(
        "--rf", "--read-format", metavar="FORMAT", dest="read_format",
        choices=READ_FORMATS, default="bed",
        help=f"Format of the read files. Support {READ_FORMATS}. Default: bed")
    parser_input.add_argument(
        "--n1", "--name1", metavar="NAME", dest="name1",
        help="Name of sample 1. If not specified, the peak file name will be "
             "used.")
    parser_input.add_argument(
        "--n2", "--name2", metavar="NAME", dest="name2",
        help="Name of sample 2. If not specified, the peak file name will be "
             "used.")

    parser_reads = parser.add_argument_group("Reads Manipulation")
    parser_reads.add_argument(
        "--s1", "--shiftsize1", metavar="N", dest="shift_size1",
        type=int, default=100,
        help="Single-end reads shift size for sample 1. Reads are shifted by "
             "`N` bp towards 3' direction and the 5' end of each shifted read "
             "is used to represent the genomic locus of the DNA fragment. "
             "Set to 0.5 * fragment size of the ChIP-seq library. "
             "Default: 100")
    parser_reads.add_argument(
        "--s2", "--shiftsize2", metavar="N", dest="shift_size2",
        type=int, default=100,
        help="Single-end reads shift size for sample 2. Default: 100")
    parser_reads.add_argument(
        "--pe", "--paired-end", dest="paired", action='store_true',
        default=False,
        help="Paired-end mode. The middle point of each read pair is used to "
             "represent the genomic locus of the DNA fragment. If specified, "
             "`--s1` and `--s2` will be ignored.")

    parser_model = parser.add_argument_group("Normalization Model Options")
    parser_model.add_argument(
        "-w", "--window-size", metavar="LENGTH", dest="window_size",
        type=_pos_int, default=2000,
        help="Window size to count reads and calculate read densities. Set to "
             "2000 is recommended for sharp histone marks like H3K4me3 or "
             "H3K27ac and 1000 for TFs or DNase-seq. Default: 2000")
    parser_model.add_argument(
        "--summit-dis", metavar="DISTANCE", dest="summit_dis_cutoff",
        type=_pos_int,
        help="Summit-to-summit distance cutoff for overlapping common peaks. "
             "Common peaks with summit distance beyond the cutoff are "
             "excluded in model fitting. Default: `window size` / 4")
    parser_model.add_argument(
        "--n-random", metavar="NUM", dest="n_random", type=int, default=10,
        help="Number of random simulations to test the enrichment of peak "
             "overlap between the specified samples. Set to 0 to disable the "
             "testing. Default: 10")

    parser_output = parser.add_argument_group("Output Options")
    parser_output.add_argument(
        "-m", "--m-cutoff", metavar="FLOAT", dest="m_cutoff",
        type=float, default=1.0,
        help="Absolute M-value (log2-ratio) cutoff to define the biased "
             "(differential binding) peaks. Default: 1.0")
    parser_output.add_argument(
        "-p", "--p-cutoff", metavar="FLOAT", dest="p_cutoff",
        type=float, default=0.01,
        help="P-value cutoff to define the biased peaks. Default: 0.01")
    parser_output.add_argument(
        "-o", "--output-dir", metavar="DIR", dest="output_dir", default=None,
        help="Output directory. Default: Current working directory")
    parser_output.add_argument(
        "--wa", "--write-all", dest="write_all", action="store_true",
        default=False,
        help="Write two extra output files containing the results of the "
             "original (unmerged) peaks.")

    parser.add_argument(
        "--verbose", dest="verbose", action="store_true", default=False,
        help="Enable verbose log messages.")
    return parser


def preprocess_args(args):
    """Pre-processing arguments."""
    args.peak_file1 = os.path.abspath(args.peak_file1)
    args.peak_file2 = os.path.abspath(args.peak_file2)
    args.read_file1 = os.path.abspath(args.read_file1)
    args.read_file2 = os.path.abspath(args.read_file2)
    args.summit_dis_cutoff = args.summit_dis_cutoff or args.window_size // 4
    args.name1 = args.name1 or os.path.splitext(
        os.path.basename(args.peak_file1))[0]
    args.name2 = args.name2 or os.path.splitext(
        os.path.basename(args.peak_file2))[0]
    args.output_dir = os.path.abspath(args.output_dir or os.getcwd())
    return args


def log_args(args):
    """Log MAnorm parameters."""
    logger.info("==== Arguments ====")
    logger.info(f"Sample 1 name = {args.name1}")
    logger.info(f"Sample 1 peak file = {args.peak_file1} [{args.peak_format}]")
    logger.info(f"Sample 1 read file = {args.read_file1} [{args.read_format}]")
    if not args.paired:
        logger.info(f"Sample 1 read shift size = {args.shift_size1}")
    logger.info(f"Sample 2 name = {args.name2}")
    logger.info(f"Sample 2 peak file = {args.peak_file2} [{args.peak_format}]")
    logger.info(f"Sample 2 read file = {args.read_file2} [{args.read_format}]")
    if not args.paired:
        logger.info(f"Sample 2 read shift size = {args.shift_size2}")
    if args.paired:
        logger.info("Paired-end mode: on")
    else:
        logger.info("Paired-end mode: off")
    logger.info(f"Window size = {args.window_size}")
    logger.info(f"Summit distance cutoff = {args.summit_dis_cutoff}")
    logger.info(f"Number of random simulation = {args.n_random}")
    logger.info(f"M-value cutoff = {args.m_cutoff}")
    logger.info(f"P-value cutoff = {args.p_cutoff}")
    logger.info(f"Output directory = {args.output_dir}")


def load_input_data(args):
    """Load required input data."""
    logger.info("Loading peaks of sample 1")
    peaks1 = load_manorm_peaks(path=args.peak_file1, format=args.peak_format,
                               name=args.name1)
    logger.info("Loading peaks of sample 2")
    peaks2 = load_manorm_peaks(path=args.peak_file2, format=args.peak_format,
                               name=args.name2)
    logger.info("Loading reads of sample 1")
    reads1 = load_reads(path=args.read_file1, format=args.read_format,
                        paired=args.paired, shift=args.shift_size1,
                        name=args.name1)
    logger.info("Loading reads of sample 2")
    reads2 = load_reads(path=args.read_file2, format=args.read_format,
                        paired=args.paired, shift=args.shift_size2,
                        name=args.name2)
    return peaks1, peaks2, reads1, reads2


def output(args, ma_model):
    """Write output files and report stats."""
    mk_dir(args.output_dir)
    if args.write_all:
        write_original_peaks(args.output_dir, ma_model.peaks1, ma_model.peaks2)
    write_all_peaks(args.output_dir, ma_model.peaks1, ma_model.peaks2,
                    ma_model.peaks_merged)
    write_wiggle_track(args.output_dir, ma_model.peaks1, ma_model.peaks2,
                       ma_model.peaks_merged)
    num_biased1, num_biased2, num_unbiased = write_biased_peaks(
        args.output_dir, ma_model.peaks1, ma_model.peaks2,
        ma_model.peaks_merged, args.m_cutoff, args.p_cutoff)
    plt_figures(args.output_dir, ma_model.peaks1, ma_model.peaks2,
                ma_model.peaks_merged, ma_model.ma_params)

    # report stats
    logger.info("==== Stats ====")
    if args.paired:
        read_type_str = 'read pairs'
    else:
        read_type_str = 'single-end reads'
    logger.info(f"Total {read_type_str} of sample 1: {ma_model.reads1.size:,}")
    logger.info(f"Total {read_type_str} of sample 2: {ma_model.reads2.size:,}")
    logger.info(
        f"Total peaks of sample 1: {ma_model.peaks1.size} "
        f"(unique: {count_unique_peaks(ma_model.peaks1)} "
        f"common: {count_common_peaks(ma_model.peaks1)})")
    logger.info(
        f"Total peaks of sample 2: {ma_model.peaks2.size} "
        f"(unique: {count_unique_peaks(ma_model.peaks2)} "
        f"common: {count_common_peaks(ma_model.peaks2)})")
    logger.info(f"Number of merged common peaks: {ma_model.peaks_merged.size}")
    logger.info(f"M-A model: M = {ma_model.ma_params[1]:.5f} * A "
                f"{ma_model.ma_params[0]:+.5f}")
    logger.info(
        f"{num_unbiased} peaks are filtered as unbiased peaks")
    logger.info(
        f"{num_biased1} peaks are filtered as sample1-biased peaks")
    logger.info(
        f"{num_biased2} peaks are filtered as sample2-biased peaks")


def run(args):
    """Run MAnorm pipeline."""
    logger.info(f"Running MAnorm {__version__}")
    log_args(args)

    logger.info("==== Running ====")
    logger.info("Step 1: Loading input data")
    peaks1, peaks2, reads1, reads2 = load_input_data(args)

    logger.info("Step 2: Processing peaks")
    ma_model = MAmodel(peaks1, peaks2, reads1, reads2)
    ma_model.process_peaks()

    logger.info("Step 3: Testing the enrichment of peak overlap")
    if args.n_random > 0:
        mean, std = random_peak_overlap(peaks1, peaks2, args.n_random)
        fc = count_common_peaks(ma_model.peaks1) / mean
        logger.info(f"Number of overlapping peaks in random: mean={mean:.2f} "
                    f"std={std:.2f}")
        logger.info(f"Fold change compared to random: {fc:.2f}")
    else:
        logger.info("Skipped")

    logger.info("Step 4: Fitting M-A normalization model on common peaks")
    ma_model.fit_model(window_size=args.window_size,
                       summit_dis_cutoff=args.summit_dis_cutoff)

    logger.info("Step 5: Normalizing all peaks")
    ma_model.normalize()

    logger.info("Step 6: Write output files")
    output(args, ma_model)


def main():
    """Main entry point, parses arguments and invoke the MAnorm application."""
    parser = configure_parser()
    args = parser.parse_args()
    args = preprocess_args(args)
    setup_logger(args.verbose)
    run(args)


if __name__ == '__main__':
    main()
