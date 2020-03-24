import pytest
from manorm.region import GenomicRegion, ManormPeak, GenomicRegions


def test_region_init():
    region = GenomicRegion(chrom='chr1', start=1, end=100)
    assert region.summit == 50
    region = GenomicRegion(chrom='chr1', start=1, end=100, summit=51)
    assert region.summit == 51
    with pytest.raises(ValueError):
        GenomicRegion(chrom='chr1', start=1, end=1)
    with pytest.raises(ValueError):
        GenomicRegion(chrom='chr1', start=2, end=1)
    with pytest.raises(ValueError):
        GenomicRegion(chrom='chr1', start=1, end=100, summit=0)
    with pytest.raises(ValueError):
        GenomicRegion(chrom='chr1', start=1, end=100, summit=100)


def test_regions_init():
    regions = GenomicRegions(name='test')
    assert regions.name == 'test'
    assert regions.size == 0
    assert regions.chroms == []


def test_regions_add():
    regions = GenomicRegions(name='test')
    region = GenomicRegion(chrom='chr1', start=1, end=100)
    regions.add(region)
    assert regions.size == 1
    assert regions.chroms == ['chr1']
    assert regions._data['chr1'] == [region]
    peak = ManormPeak(chrom='chr1', start=1, end=100)
    regions.add(peak)
    assert regions.size == 2
    assert regions.chroms == ['chr1']
    assert regions._data['chr1'] == [region, peak]
    with pytest.raises(ValueError):
        regions.add(123)


def test_regions_fetch():
    regions = GenomicRegions(name='test')
    region = GenomicRegion(chrom='chr1', start=1, end=100)
    regions.add(region)
    assert regions.size == 1
    assert regions.chroms == ['chr1']
    assert regions.fetch(chrom='chr1') == [region]
    assert regions.fetch(chrom='chr11') == []


def test_regions_sort():
    regions = GenomicRegions(name='test')
    region1 = GenomicRegion(chrom='chr1', start=1, end=100)
    region2 = GenomicRegion(chrom='chr1', start=10, end=50)
    regions.add(region=region1)
    regions.add(region=region2)
    assert regions.size == 2
    assert regions.chroms == ['chr1']
    regions.sort(by='start', ascending=True)
    assert regions.fetch(chrom='chr1')[0] == region1
    assert regions.fetch(chrom='chr1')[1] == region2
    regions.sort(by='start', ascending=False)
    assert regions.fetch(chrom='chr1')[0] == region2
    assert regions.fetch(chrom='chr1')[1] == region1
    regions.sort(by='end', ascending=True)
    assert regions.fetch(chrom='chr1')[0] == region2
    assert regions.fetch(chrom='chr1')[1] == region1
    regions.sort(by='end', ascending=False)
    assert regions.fetch(chrom='chr1')[0] == region1
    assert regions.fetch(chrom='chr1')[1] == region2
    regions.sort(by='summit', ascending=True)
    assert regions.fetch(chrom='chr1')[0] == region2
    assert regions.fetch(chrom='chr1')[1] == region1
    regions.sort(by='summit', ascending=False)
    assert regions.fetch(chrom='chr1')[0] == region1
    assert regions.fetch(chrom='chr1')[1] == region2
