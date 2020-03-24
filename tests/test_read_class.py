import pytest

from manorm.read import Reads


def test_reads_init():
    reads = Reads(name='test')
    assert reads.name == 'test'
    assert reads.size == 0
    assert reads.chroms == []


def test_reads_add():
    reads = Reads(name='test')
    reads.add('chr1', 100)
    assert reads.size == 1
    assert reads.chroms == ['chr1']
    assert reads._data['chr1'] == [100]


def test_reads_sort():
    reads = Reads(name='test')
    reads.add('chr1', 100)
    reads.add('chr1', 102)
    reads.add('chr1', 1)
    reads.sort()
    assert reads.size == 3
    assert reads.chroms == ['chr1']
    assert reads._data['chr1'] == [1, 100, 102]


def test_reads_count():
    reads = Reads(name='test')
    reads.add('chr1', 100)
    reads.add('chr1', 102)
    reads.add('chr1', 1)
    reads.sort()
    with pytest.raises(ValueError):
        reads.count('chr1', 1, 1)
    assert reads.count('chr11', 1, 200) == 0
    assert reads.count('chr1', -100, 0) == 0
    assert reads.count('chr1', 1, 2) == 1
    assert reads.count('chr1', 1, 100) == 1
    assert reads.count('chr1', 1, 101) == 2
    assert reads.count('chr1', 1, 200) == 3
