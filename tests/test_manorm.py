import os

import pytest

from manorm import application


def test_manorm(data_dir, tmp_dir):
    peak_file1 = os.path.join(data_dir, 'H1hescH3k4me3Rep1_peaks.xls')
    peak_file2 = os.path.join(data_dir, 'K562H3k4me3Rep1_peaks.xls')
    read_file1 = os.path.join(data_dir, 'H1hescH3k4me3Rep1_reads.bed')
    read_file2 = os.path.join(data_dir, 'K562H3k4me3Rep1_reads.bed')
    shift_size1 = 100
    shift_size2 = 100
    window_size = 2000
    summit_dis = 500
    n_random = 5
    m_cutoff = 1
    p_cutoff = 0.01
    app = application.MAnorm(peak_file1=peak_file1, peak_file2=peak_file2, read_file1=read_file1, read_file2=read_file2,
                             peak_format='macs', read_format='bed', name1='H1_H3K4me3', name2='K562_H3K4me3',
                             shift_size1=shift_size1, shift_size2=shift_size2, paired=False,
                             window_size=window_size, summit_dis=summit_dis, n_random=n_random,
                             m_cutoff=m_cutoff, p_cutoff=p_cutoff, write_all=True,
                             output_dir=tmp_dir)
    app.run()
    fn1 = os.path.join(data_dir, 'H1_H3K4me3_vs_K562_H3K4me3_all_MAvalues.xls')
    fn2 = os.path.join(tmp_dir, 'H1_H3K4me3_vs_K562_H3K4me3_all_MAvalues.xls')
    with open(fn1, "r") as f1, open(fn2, "r") as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
        assert len(lines1) == len(lines2)
        for idx, tmp_line1 in enumerate(lines1):
            tmp_line1 = tmp_line1.strip().split('\t')
            tmp_line2 = lines2[idx].strip().split('\t')
            if idx == 0:
                assert tmp_line1 == tmp_line2
            else:
                assert tmp_line1[:4] == tmp_line2[:4]
                assert float(tmp_line1[4]) == pytest.approx(float(tmp_line2[4]), abs=1e-5)
                assert float(tmp_line1[5]) == pytest.approx(float(tmp_line2[5]), abs=1e-5)
                assert float(tmp_line1[6]) == pytest.approx(float(tmp_line2[6]), abs=1e-5)
                assert tmp_line1[7] == tmp_line2[7]
                assert float(tmp_line1[8]) == pytest.approx(float(tmp_line2[8]), abs=1e-5)
                assert float(tmp_line1[9]) == pytest.approx(float(tmp_line2[9]), abs=1e-5)
