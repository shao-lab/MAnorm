import os

import pytest

from manorm.read import load_reads


@pytest.fixture
def data_dir(tmpdir):
    """Return the directory of data files."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


def test_bed(data_dir):
    reads = load_reads(os.path.join(data_dir, 'test_reads.bed'), format='bed', paired=False, shift=100)
    assert sorted(reads.chroms) == ['chr1', 'chr2', 'chr9']
    assert reads.size == 4
    assert reads.data['chr1'] == [101, 400]
    assert reads.data['chr2'] == [12445]
    assert reads.data['chr9'] == [12245]


def test_bedpe(data_dir):
    reads = load_reads(os.path.join(data_dir, 'test_reads.bedpe'), format='bedpe', paired=True)
    assert sorted(reads.chroms) == ['chr1', 'chr9']
    assert reads.size == 2
    assert reads.data['chr1'] == [200]
    assert reads.data['chr9'] == [2400]


def test_sam(data_dir):
    reads = load_reads(os.path.join(data_dir, 'test_reads.sam'), format='sam', paired=False, shift=100)
    assert sorted(reads.chroms) == ['chr1', 'chr2', 'chr9']
    assert reads.size == 4
    assert reads.data['chr1'] == [101, 400]
    assert reads.data['chr2'] == [12445]
    assert reads.data['chr9'] == [12245]


def test_bam(data_dir):
    reads = load_reads(os.path.join(data_dir, 'test_reads.bam'), format='bam', paired=False, shift=100)
    assert sorted(reads.chroms) == ['chr1', 'chr2', 'chr9']
    assert reads.size == 4
    assert reads.data['chr1'] == [101, 400]
    assert reads.data['chr2'] == [12445]
    assert reads.data['chr9'] == [12245]


def test_bam_pe(data_dir):
    reads = load_reads(os.path.join(data_dir, 'test_reads_pe.bam'), format='bam', paired=True)
    assert sorted(reads.chroms) == ['chr1', 'chr9']
    assert reads.size == 2
    assert reads.data['chr1'] == [150]
    assert reads.data['chr9'] == [2000]
