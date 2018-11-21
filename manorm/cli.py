# -*- coding: utf-8 -*-

"""
manorm.cli
~~~~~~~~~~

Command line implementation of MAnorm.
"""

from __future__ import absolute_import

import argparse
import os
from textwrap import dedent

from manorm import __version__, application
from manorm.peak import PEAK_FORMATS
from manorm.read import READ_FORMATS


def argparser_config():
    """Configure the arguments parser."""
    description = """MAnorm -- A robust model for quantitative comparison of ChIP-seq data sets."""
    epilog = dedent("""
    Citation:
    
    Shao Z, Zhang Y, Yuan GC, Orkin SH, Waxman DJ. MAnorm: a robust model for quantitative 
    comparison of ChIP-Seq data sets. Genome biology. 2012 Mar 16;13(3):R16.
    https://doi.org/10.1186/gb-2012-13-3-r16
    
    For more information, please visit:
    Documentation: https://manorm.readthedocs.io
    Source code: https://github.com/shao-lab/MAnorm
    Report bugs: https://github.com/shao-lab/MAnorm/issues
    """)

    parser = argparse.ArgumentParser(description=description, epilog=epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--version", action="version", version="MAnorm {}".format(__version__))

    input_args = parser.add_argument_group("Input Arguments")
    input_args.add_argument("--p1", "--peak1", metavar="FILE", dest="peak_file1", type=existed_file, required=True,
                            help="Peak file of sample 1.")
    input_args.add_argument("--p2", "--peak2", metavar="FILE", dest="peak_file2", type=existed_file, required=True,
                            help="Peak file of sample 2.")
    input_args.add_argument("--pf", "--peak-format", metavar="FORMAT", dest="peak_format",
                            choices=PEAK_FORMATS, default="bed",
                            help="Format of peak files. Supported: {}. Default: bed".format(PEAK_FORMATS))
    input_args.add_argument("--r1", "--read1", metavar="FILE", dest="read_file1", type=existed_file, required=True,
                            help="Read file of sample 1.")
    input_args.add_argument("--r2", "--read2", metavar="FILE", dest="read_file2", type=existed_file, required=True,
                            help="Read file of sample 2.")
    input_args.add_argument("--rf", "--read-format", metavar="FORMAT", dest="read_format",
                            choices=READ_FORMATS, default="bed",
                            help="Format of read files. Supported: {}. Default: bed".format(READ_FORMATS))
    input_args.add_argument("--n1", "--name1", metavar="NAME", dest="name1", default=None,
                            help="Name of sample 1 (i.e. experiment condition/cell-type), which is used in output "
                                 "files. If not specified, peak/read file names will be used as sample names.")
    input_args.add_argument("--n2", "--name2", metavar="NAME", dest="name2", default=None,
                            help="Name of sample 2.")

    reads_args = parser.add_argument_group("Reads Manipulation")
    reads_args.add_argument("--s1", "--shiftsize1", metavar="N", dest="shift_size1", type=int, default=100,
                            help="Reads shift size of sample 1. Shift single-end reads by N bp towards 3' direction "
                                 "and the 5' end of each shifted read is used to represent the genomic locus of "
                                 "underlying DNA fragment. Set to half of DNA fragment size of the ChIP-seq library. "
                                 "Default: 100")
    reads_args.add_argument("--s2", "--shiftsize2", metavar="N", dest="shift_size2", type=int, default=100,
                            help="Reads shift size of sample 2. Default: 100")
    reads_args.add_argument("--paired", dest="paired", action='store_true', default=False,
                            help="Paired-end mode. The middle point of each read pair is used to represent the genomic "
                                 "locus of underlying DNA fragment. `--s1`, `--s2` are ignored with this option on.")

    model_args = parser.add_argument_group("Normalization Model")
    model_args.add_argument("-w", "--window-size", metavar="LENGTH", dest="window_size", type=pos_int, default=2000,
                            help="Window size to count reads and calculate read densities. 2000 is recommended for "
                                 "sharp histone marks like H3K4me3 and H3K27ac, and 1000 for TFs or DNase-seq. "
                                 "Default: 2000")
    model_args.add_argument("--summit-dis", metavar="LENGTH", dest="summit_dis", type=pos_int, default=None,
                            help="Overlapping common peaks with summit-to-summit distance above this are excluded "
                                 "in model fitting. Default: `-w/--window-size` / 4")
    model_args.add_argument("--n-random", metavar="N", dest="n_random", type=pos_int, default=10,
                            help="Number of random simulations to test the enrichment of peak overlap between two "
                                 "samples. Set to 0 to disable the testing. Default: 10")

    output_args = parser.add_argument_group("Output Arguments")
    output_args.add_argument("-m", "--m-cutoff", metavar="FLOAT", dest="m_cutoff", type=float, default=1.0,
                             help="Absolute M-value (log2-ratio) cutoff to define biased (differential binding) peaks. "
                                  "Default: 1.0")
    output_args.add_argument("-p", "--p-cutoff", metavar="FLOAT", dest="p_cutoff", type=float, default=0.01,
                             help="P-value cutoff to define biased peaks. Default: 0.01")
    output_args.add_argument("--write-all", dest="write_all", action="store_true", default=False,
                             help="Write two extra output files containing the results of original (unmerged) peaks.")
    output_args.add_argument("--verbose", dest="verbose", action="store_true", default=False,
                             help="Enable verbose log messages.")
    output_args.add_argument("--dir", "--output-dir", metavar="DIR", dest="output_dir", default=None,
                             help="Output directory path. Default: Current working directory")
    return parser


def existed_file(path):
    """Wrapper function to check whether a passed argument is an existed file."""
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError("file not found: {!r}".format(path))
    return path


def pos_int(value):
    """wrapper function to check whether a passed argument is a positive integer."""
    try:
        value_int = int(value)
        if value_int <= 0:
            raise ValueError
    except ValueError:
        raise argparse.ArgumentTypeError("invalid positive int value: {!r}".format(value))
    return value_int


def main():
    """MAnorm main entry point.
    This function parses arguments and passes them to an instance of :class:`MAnorm` to run MAnorm application.
    """
    parser = argparser_config()
    args = parser.parse_args()
    app = application.MAnorm(args)
    app.run()


if __name__ == '__main__':
    main()
