"""MAnorm main script for running from the command line."""

import argparse
from manorm import __version__, workflow


def argparser_config():
    """Configure the arguments parser.
    """
    description = """MAnorm -- A robust model for quantitative comparison of ChIP-Seq data sets."""

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", action="version", version="MAnorm {}".format(__version__))

    input_args = parser.add_argument_group("Input Arguments")
    input_args.add_argument("--p1", dest="peaks_file1", type=str, required=True,
                            help="Peaks file of sample 1. BED and MACS format are currently supported. Please refer to "
                                 "documents for details.")
    input_args.add_argument("--p2", dest="peaks_file2", type=str, required=True,
                            help="Peaks file of sample 2.")
    input_args.add_argument("--r1", dest="reads_file1", type=str, required=True,
                            help="Reads file of sample 1. BED format are currently supported.")
    input_args.add_argument("--r2", dest="reads_file2", type=str, required=True,
                            help="Reads file of sample 2.")
    input_args.add_argument("--s1", dest="shift_size1", type=int, default=100,
                            help="Reads shift size of sample 1. This value is used to shift reads towards 3' direction "
                                 "to determine the precise binding site. Set as half of the fragment length.")
    input_args.add_argument("--s2", dest="shift_size2", type=int, default=100,
                            help="Reads shift size of sample 2.")

    model_args = parser.add_argument_group("Normalization Model Arguments")
    model_args.add_argument("-w", dest="width", type=int, default=1000,
                            help="Width of the window to calculate read densities. Windows with unified length of "
                                 "2*width centered at peak summit/midpoint are used to qualify the binding signal. "
                                 "This should match the typical length of peaks, a value of 1000 is recommended for "
                                 "sharp histone marks like H3K4me3 and H3K9/27ac, and 500 for transcription factors "
                                 "or DNase-Seq.")
    model_args.add_argument("-d", dest="dis_cutoff", type=int, default=None,
                            help="Summit-to-summit distance cutoff for common peaks. Default=width/2. Only overlapped "
                                 "peaks with summit-to-summit distance less than than this value are considered as "
                                 "real common peaks of two samples when fitting M-A normalization model.")

    advanced_args = parser.add_argument_group("Advanced Arguments")
    advanced_args.add_argument("-n", dest="n_random", type=int, default=10,
                               help="Number of simulation to test the enrichment of peak overlap between two samples.")
    advanced_args.add_argument("-m", dest="m_cutoff", type=float, default=1.0,
                               help="M-value cutoff to distinguish biased peaks from unbiased peaks. Peaks with "
                                    "M-value>=M_cutoff and P-value<=P_cutoff are defined as sample1-biased(specific) "
                                    "peaks, while peaks with M-value<=-1*M_cutoff and P-value<=P_cutoff are defined "
                                    "as sample2-biased peaks.")
    advanced_args.add_argument("-p", dest="p_cutoff", type=float, default=0.01,
                               help="P-value cutoff to define biased (sample 1/2-specific) peaks.")

    output_args = parser.add_argument_group("Output arguments")
    output_args.add_argument("-s", dest="full_output", action="store_true", default=False,
                             help="By default, MAnorm will write the comparison results of unique and merged common "
                                  "peaks in a single output file. With this option on, two extra files which contains "
                                  "the results of the original(unmerged) peaks will also be output.")
    output_args.add_argument("--name1", dest="name1", type=str, default=None,
                             help="Name (experiment condition/cell-type etc.) of sample1. If specified, it will be "
                                  "used to replace the peaks/reads input file name as the sample name in output files.")
    output_args.add_argument("--name2", dest="name2", type=str, default=None,
                             help="Name (experiment condition/cell-type etc.) of sample2.")
    output_args.add_argument("-o", dest="output_name", type=str, required=True,
                             help="Output directory. When --name1 and --name2 are not specified, MAnorm will use it as "
                                  "the prefix of comparison output file.")
    return parser


def main():
    """Main entry point for MAnorm."""
    parser = argparser_config()
    args = parser.parse_args()
    peaks_file1 = args.peaks_file1
    peaks_file2 = args.peaks_file2
    reads_file1 = args.reads_file1
    reads_file2 = args.reads_file2
    shiftsize1 = args.shift_size1
    shiftsize2 = args.shift_size2
    peak_width = args.width
    if args.dis_cutoff:
        summit_dis_cutoff = args.dis_cutoff
    else:
        summit_dis_cutoff = peak_width / 2
    n_random = args.n_random
    m_cutoff = args.m_cutoff
    p_cutoff = args.p_cutoff
    full_output = args.full_output
    name1 = args.name1
    name2 = args.name2
    output_name = args.output_name
    workflow.main(peaks_file1=peaks_file1, peaks_file2=peaks_file2, reads_file1=reads_file1, reads_file2=reads_file2,
                  shift_size1=shiftsize1, shift_size2=shiftsize2, peak_width=peak_width,
                  summit_dis_cutoff=summit_dis_cutoff, n_random=n_random, m_cutoff=m_cutoff, p_cutoff=p_cutoff,
                  full_output=full_output, name1=name1, name2=name2, output_name=output_name)


if __name__ == '__main__':
    main()
