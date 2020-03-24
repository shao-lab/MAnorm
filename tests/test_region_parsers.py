import os

from manorm.region import load_manorm_peaks


def test_bed_parser(data_dir):
    peaks = load_manorm_peaks(os.path.join(data_dir, 'test_peaks.bed'),
                              format='bed')
    assert peaks.chroms == ['chr1', 'chr2', 'chr9']
    assert peaks.size == 4
    assert peaks.fetch('chr1')[0].start == 1
    assert peaks.fetch('chr1')[0].end == 100
    assert peaks.fetch('chr1')[0].summit == 50
    assert peaks.fetch('chr1')[1].start == 2
    assert peaks.fetch('chr1')[1].end == 200
    assert peaks.fetch('chr1')[1].summit == 101
    assert peaks.fetch('chr2')[0].start == 1
    assert peaks.fetch('chr2')[0].end == 150
    assert peaks.fetch('chr2')[0].summit == 75
    assert peaks.fetch('chr9')[0].start == 5
    assert peaks.fetch('chr9')[0].end == 123
    assert peaks.fetch('chr9')[0].summit == 64


def test_bed3_summit_parser(data_dir):
    peaks = load_manorm_peaks(os.path.join(data_dir, 'test_peaks_summit.bed'),
                              format='bed3-summit')
    assert sorted(peaks.chroms) == ['chr1', 'chr2', 'chr9']
    assert peaks.size == 4
    assert peaks.fetch('chr1')[0].start == 1
    assert peaks.fetch('chr1')[0].end == 100
    assert peaks.fetch('chr1')[0].summit == 50
    assert peaks.fetch('chr1')[1].start == 2
    assert peaks.fetch('chr1')[1].end == 200
    assert peaks.fetch('chr1')[1].summit == 100
    assert peaks.fetch('chr2')[0].start == 1
    assert peaks.fetch('chr2')[0].end == 150
    assert peaks.fetch('chr2')[0].summit == 2
    assert peaks.fetch('chr9')[0].start == 5
    assert peaks.fetch('chr9')[0].end == 123
    assert peaks.fetch('chr9')[0].summit == 55


def test_macs_parser(data_dir):
    peaks = load_manorm_peaks(os.path.join(data_dir, 'test_peaks_macs.xls'),
                              format='macs')
    assert sorted(peaks.chroms) == ['chr1', 'chr2', 'chr22']
    assert peaks.size == 9
    assert peaks.fetch('chr1')[0].start == 16192292
    assert peaks.fetch('chr1')[0].end == 16193176
    assert peaks.fetch('chr1')[0].summit == 16192491
    assert peaks.fetch('chr1')[1].start == 17081409
    assert peaks.fetch('chr1')[1].end == 17082059
    assert peaks.fetch('chr1')[1].summit == 17081819
    assert peaks.fetch('chr2')[0].start == 17082916
    assert peaks.fetch('chr2')[0].end == 17084523
    assert peaks.fetch('chr2')[0].summit == 17084177
    assert peaks.fetch('chr22')[0].start == 17565233
    assert peaks.fetch('chr22')[0].end == 17567384
    assert peaks.fetch('chr22')[0].summit == 17565935


def test_macs2_parser(data_dir):
    peaks = load_manorm_peaks(os.path.join(data_dir, 'test_peaks_macs2.xls'),
                              format='macs2')
    assert sorted(peaks.chroms) == ['chr1', 'chr2', 'chr22']
    assert peaks.size == 10
    assert peaks.fetch('chr1')[0].start == 569795
    assert peaks.fetch('chr1')[0].end == 570052
    assert peaks.fetch('chr1')[0].summit == 569927
    assert peaks.fetch('chr1')[1].start == 713873
    assert peaks.fetch('chr1')[1].end == 714348
    assert peaks.fetch('chr1')[1].summit == 714069
    assert peaks.fetch('chr2')[0].start == 778179
    assert peaks.fetch('chr2')[0].end == 778484
    assert peaks.fetch('chr2')[0].summit == 778368
    assert peaks.fetch('chr22')[0].start == 834127
    assert peaks.fetch('chr22')[0].end == 834359
    assert peaks.fetch('chr22')[0].summit == 834280


def test_narrowpeak_parser(data_dir):
    peaks = load_manorm_peaks(os.path.join(data_dir, 'test_peaks.narrowPeak'),
                              format='narrowpeak')
    assert sorted(peaks.chroms) == ['chr1', 'chr2', 'chr22']
    assert peaks.size == 10
    assert peaks.fetch('chr1')[0].start == 569795
    assert peaks.fetch('chr1')[0].end == 570052
    assert peaks.fetch('chr1')[0].summit == 569927
    assert peaks.fetch('chr1')[1].start == 713873
    assert peaks.fetch('chr1')[1].end == 714348
    assert peaks.fetch('chr1')[1].summit == 714069
    assert peaks.fetch('chr2')[0].start == 778179
    assert peaks.fetch('chr2')[0].end == 778484
    assert peaks.fetch('chr2')[0].summit == 778368
    assert peaks.fetch('chr22')[0].start == 834127
    assert peaks.fetch('chr22')[0].end == 834359
    assert peaks.fetch('chr22')[0].summit == 834280
