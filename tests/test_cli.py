import os
from argparse import Namespace

import pytest

from manorm.cli import configure_parser, preprocess_args, run


def test_configure_parser(data_dir):
    parser = configure_parser()
    args = parser.parse_args(
        ["--p1", os.path.join(data_dir, 'H1hescH3k4me3Rep1_peaks.xls'),
         "--p2", os.path.join(data_dir, 'K562H3k4me3Rep1_peaks.xls'),
         "--r1", os.path.join(data_dir, 'H1hescH3k4me3Rep1_reads.bed'),
         "--r2", os.path.join(data_dir, 'K562H3k4me3Rep1_reads.bed')])
    args = preprocess_args(args)
    assert args.name1 == "H1hescH3k4me3Rep1_peaks"
    assert args.name2 == "K562H3k4me3Rep1_peaks"
    assert args.summit_dis_cutoff == args.window_size // 4
    assert args.output_dir == os.getcwd()


def test_run(data_dir, tmp_dir):
    args = Namespace(**dict(
        peak_file1=os.path.join(data_dir, 'H1hescH3k4me3Rep1_peaks.xls'),
        peak_file2=os.path.join(data_dir, 'K562H3k4me3Rep1_peaks.xls'),
        read_file1=os.path.join(data_dir, 'H1hescH3k4me3Rep1_reads.bed'),
        read_file2=os.path.join(data_dir, 'K562H3k4me3Rep1_reads.bed'),
        peak_format='macs', read_format='bed',
        name1='H1_H3K4me3', name2='K562_H3K4me3',
        shift_size1=100, shift_size2=100, paired=False,
        window_size=2000, summit_dis_cutoff=500, n_random=5,
        m_cutoff=1, p_cutoff=0.01, write_all=True, output_dir=tmp_dir))
    run(args)
    fn1 = os.path.join(data_dir, 'H1_H3K4me3_vs_K562_H3K4me3_all_MAvalues.xls')
    fn2 = os.path.join(tmp_dir, 'H1_H3K4me3_vs_K562_H3K4me3_all_MAvalues.xls')
    with open(fn1, "r") as f1, open(fn2, "r") as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
        assert len(lines1) == len(lines2)
        assert lines1[0] == lines2[0]
        for tmp_line1, tmp_line2 in zip(lines1[1:], lines2[1:]):
            tmp_line1 = tmp_line1.strip().split('\t')
            tmp_line2 = tmp_line2.strip().split('\t')
            assert tmp_line1[:4] == tmp_line2[:4]
            assert float(tmp_line1[4]) == pytest.approx(
                float(tmp_line2[4]), abs=1e-5)
            assert float(tmp_line1[5]) == pytest.approx(
                float(tmp_line2[5]), abs=1e-5)
            assert float(tmp_line1[6]) == pytest.approx(
                float(tmp_line2[6]), abs=1e-5)
            assert tmp_line1[7] == tmp_line2[7]
            assert float(tmp_line1[8]) == pytest.approx(
                float(tmp_line2[8]), abs=1e-5)
            assert float(tmp_line1[9]) == pytest.approx(
                float(tmp_line2[9]), abs=1e-5)
