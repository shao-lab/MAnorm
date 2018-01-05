import os
import shutil
import manorm
from manorm import workflow


def setup():
    data_dir = os.path.join(os.path.dirname(os.path.abspath(manorm.__file__)), 'tests', 'data')
    os.chdir(data_dir)


def teardown():
    shutil.rmtree('H1_vs_K562_H3K4me3')


def test_manorm():
    peaks_file1 = 'H1hescH3k4me3Rep1_peaks.xls'
    peaks_file2 = 'K562H3k4me3Rep1_peaks.xls'
    reads_file1 = 'H1hescH3k4me3Rep1_reads.bed'
    reads_file2 = 'K562H3k4me3Rep1_reads.bed'
    shiftsize1 = 100
    shiftsize2 = 100
    peak_width = 1000
    summit_dis_cutoff = 500
    n_random = 5
    m_cutoff = 1
    p_cutoff = 0.01
    full_output = False
    name1 = None
    name2 = None
    output_name = 'H1_vs_K562_H3K4me3'
    workflow.main(peaks_file1=peaks_file1, peaks_file2=peaks_file2, reads_file1=reads_file1, reads_file2=reads_file2,
                  shift_size1=shiftsize1, shift_size2=shiftsize2, peak_width=peak_width,
                  summit_dis_cutoff=summit_dis_cutoff, n_random=n_random, m_cutoff=m_cutoff, p_cutoff=p_cutoff,
                  full_output=full_output, name1=name1, name2=name2, output_name=output_name)

    fn1 = os.path.join('example', 'H1_vs_K562_H3K4me3_all_MAvalues.xls')
    fn2 = os.path.join('H1_vs_K562_H3K4me3', 'H1_vs_K562_H3K4me3_all_MAvalues.xls')
    with open(fn1, "r") as f1, open(fn2, "r") as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
        assert lines1 == lines2
