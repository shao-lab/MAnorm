import os

from manorm.read import load_reads


def test_bed(data_dir):
    reads = load_reads(os.path.join(data_dir, 'test_reads.bed'), format='bed',
                       paired=False, shift=100)
    assert reads.chroms == ['chr1', 'chr2', 'chr9']
    assert reads.size == 4
    assert reads._data['chr1'] == [101, 400]
    assert reads._data['chr2'] == [12445]
    assert reads._data['chr9'] == [12245]


def test_bedpe(data_dir):
    reads = load_reads(os.path.join(data_dir, 'test_reads.bedpe'),
                       format='bedpe', paired=True)
    assert reads.chroms == ['chr1', 'chr9']
    assert reads.size == 2
    assert reads._data['chr1'] == [200]
    assert reads._data['chr9'] == [2400]


def test_sam(data_dir):
    reads = load_reads(os.path.join(data_dir, 'test_reads.sam'), format='sam',
                       paired=False, shift=100)
    assert reads.chroms == ['chr1', 'chr2', 'chr9']
    assert reads.size == 4
    assert reads._data['chr1'] == [101, 400]
    assert reads._data['chr2'] == [12445]
    assert reads._data['chr9'] == [12245]


def test_bam(data_dir):
    reads = load_reads(os.path.join(data_dir, 'test_reads.bam'), format='bam',
                       paired=False, shift=100)
    assert reads.chroms == ['chr1', 'chr2', 'chr9']
    assert reads.size == 4
    assert reads._data['chr1'] == [101, 400]
    assert reads._data['chr2'] == [12445]
    assert reads._data['chr9'] == [12245]


def test_bam_pe(data_dir):
    reads = load_reads(os.path.join(data_dir, 'test_reads_pe.bam'),
                       format='bam', paired=True)
    assert reads.chroms == ['chr1', 'chr9']
    assert reads.size == 2
    assert reads._data['chr1'] == [150]
    assert reads._data['chr9'] == [2000]
