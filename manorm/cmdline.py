"""MAnorm main script for running from the command line."""
import os
import argparse
from manorm import workflow


def argparser_config():
    """Configure the arguments parser.
    """
    description = """MAnorm -- A robust model for quantitative comparison of ChIP-Seq data sets."""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    group_input = parser.add_argument_group("Input File Arguments")
    group_input.add_argument("--p1", dest="peaks_file1", type=str, required=True,
                             help="Path of peaks file of sample 1. BED and MACS format are currently supported."
                                  "Please refer to documents for details.")
    group_input.add_argument("--p2", dest="peaks_file2", type=str, required=True,
                             help="Path of peaks file of sample 2.")
    group_input.add_argument("--r1", dest="reads_file1", type=str, required=True,
                             help="Path of reads file of sample 1. BED format are currently supported.")
    group_input.add_argument("--r2", dest="reads_file2", type=str, required=True,
                             help="Path of reads file of sample 2.")
    group_input.add_argument("--s1", dest="shift_size1", type=int, default=100,
                             help="Reads shiftsize of sample 1. This value is used to shift reads towards 3' direction"
                                  "to account for the actual binding site. Set as half of the fragment length.")
    group_input.add_argument("--s2", dest="shift_size2", type=int, default=100,
                             help="Reads shiftsize of sample 2.")

    group_model = parser.add_argument_group("Model arguments")
    group_model.add_argument("-w", dest="width", type=int, default=1000,
                             help="Half width of the window size when calculating read densities. Each window with "
                                  "length of 2*width is centered at peak summit or midpoint. This value should match "
                                  "the typical length of peaks, thus we recommend 1000 for sharp histone marks like "
                                  "H3K4me3 and H3K9/27ac, or 500 for transcription factors or DNase-Seq.")
    group_model.add_argument("-d", dest="distance_cutoff", type=int,
                             help="Summit to summit distance cutoff for common peaks. Default=width/2. Only peaks "
                                  "with summit distance less than than this value are considered as real common peaks "
                                  "of two samples.")

    group_advanced = parser.add_argument_group("Advanced arguments")
    group_advanced.add_argument("-n", dest="random_times", type=int, default=5,
                                help="Times of permutation to test the enrichment of peak overlap between two samples.")
    group_advanced.add_argument("-v", dest="overlap_dependent", action="store_true", default=False,
                                help="With this option on, MAnorm will pick out biased/unbiased peaks in an "
                                     "overlap-dependent manner. Biased peaks are only chosen from unique peaks, "
                                     "and unbiased peaks are only chosen from common peaks.")
    group_advanced.add_argument("-p", dest="p_cutoff", type=float, default=0.01,
                                help="P-value cutoff to define biased (sample 1/2-specific) peaks.")
    group_advanced.add_argument("-m", dest="m_cutoff", type=float, default=1.0,
                                help="M-value cutoff to distinguish biased peaks from unbiased peaks. Peaks with "
                                     "M-value>=M_cutoff and P-value<=P_cutoff as defined as sample1-biased(specific) "
                                     "peaks, while peaks with M-value<=-1*M_cutoff and P-value<=P_cutoff are defined "
                                     "as sample2-biased peaks.")

    group_output = parser.add_argument_group("Output arguments")
    group_output.add_argument("-s", dest="output_all", action="store_true", default=False,
                              help="By default, MAnorm will output the results of unique and merged common peaks of "
                                   "two samples. With this option on, MAnorm will output two extra files containing"
                                   "the results of the original(unmerged) peaks in two samples.")
    group_output.add_argument("-o", dest="output_prefix", type=str, required=True,
                              help="Comparison name, this is used as the folder name and prefix of output files.")
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
    if args.distance_cutoff:
        distance_cutoff = args.distance_cutoff
    else:
        distance_cutoff = peak_width / 2
    random_times = args.random_times
    overlap_dependent = args.overlap_dependent
    p_cutoff = args.p_cutoff
    m_cutoff = args.m_cutoff
    output_all = args.output_all
    output_name = args.output_prefix

    workflow.main(peaks_file1=peaks_file1, peaks_file2=peaks_file2, reads_file1=reads_file1, reads_file2=reads_file2,
                  shift_size1=shiftsize1, shift_size2=shiftsize2, peak_width=peak_width, distance_cutoff=distance_cutoff,
                  random_times=random_times, overlap_dependent=overlap_dependent, m_cutoff=m_cutoff, p_cutoff=p_cutoff,
                  output_all=output_all, output_name=output_name)


if __name__ == '__main__':
    main()
