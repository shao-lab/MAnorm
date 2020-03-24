from manorm.region import GenomicRegion, ManormPeak, GenomicRegions
from manorm.region.utils import overlap_on_same_chrom, \
    classify_peaks_by_overlap, merge_common_peaks, generate_random_regions, \
    count_common_peaks, count_unique_peaks


def test_region_overlap_on_same_chrom():
    peak1 = GenomicRegion(chrom='chr1', start=1, end=100, summit=50)
    peak2 = GenomicRegion(chrom='chr1', start=50, end=150, summit=120)
    peak3 = GenomicRegion(chrom='chr1', start=150, end=200, summit=166)
    peak4 = GenomicRegion(chrom='chr1', start=250, end=500, summit=350)
    flag1, flag2 = overlap_on_same_chrom([peak1, peak3], [peak2, peak4])
    assert flag1[0]
    assert not flag1[1]
    assert flag2[0]
    assert not flag2[1]


def test_classify_peak_overlap():
    peaks1 = GenomicRegions(name='test1')
    peaks2 = GenomicRegions(name='test2')
    peak1 = ManormPeak(chrom='chr1', start=1, end=100, summit=50)
    peak2 = ManormPeak(chrom='chr1', start=50, end=150, summit=120)
    peak3 = ManormPeak(chrom='chr1', start=150, end=200, summit=166)
    peak4 = ManormPeak(chrom='chr1', start=250, end=500, summit=350)
    peak5 = ManormPeak(chrom='chr2', start=300, end=400, summit=333)
    peak6 = ManormPeak(chrom='chr2', start=1, end=100, summit=50)
    peak7 = ManormPeak(chrom='chr22', start=1, end=100, summit=50)
    peaks1.add(peak1)
    peaks1.add(peak3)
    peaks1.add(peak5)
    peaks1.add(peak7)
    peaks2.add(peak2)
    peaks2.add(peak4)
    peaks2.add(peak6)
    peaks1, peaks2 = classify_peaks_by_overlap(peaks1, peaks2)
    assert peaks1.fetch('chr1')[0].iscommon
    assert not peaks1.fetch('chr1')[1].iscommon
    assert peaks2.fetch('chr1')[0].iscommon
    assert not peaks2.fetch('chr1')[1].iscommon
    assert not peaks1.fetch('chr2')[0].iscommon
    assert not peaks2.fetch('chr2')[0].iscommon
    assert not peaks1.fetch('chr22')[0].iscommon


def test_merge_common_peaks():
    peaks1 = GenomicRegions(name='test1')
    peaks2 = GenomicRegions(name='test2')
    peak1 = ManormPeak(chrom='chr1', start=1, end=100, summit=50)
    peak2 = ManormPeak(chrom='chr1', start=50, end=150, summit=120)
    peak3 = ManormPeak(chrom='chr1', start=149, end=200, summit=166)
    peak4 = ManormPeak(chrom='chr1', start=250, end=500, summit=350)
    peak5 = ManormPeak(chrom='chr2', start=300, end=400, summit=333)
    peak6 = ManormPeak(chrom='chr2', start=1, end=100, summit=50)
    peak7 = ManormPeak(chrom='chr2', start=50, end=120, summit=80)
    peak8 = ManormPeak(chrom='chr22', start=1, end=100, summit=50)
    peaks1.add(peak1)
    peaks1.add(peak3)
    peaks1.add(peak5)
    peaks1.add(peak7)
    peaks2.add(peak2)
    peaks2.add(peak4)
    peaks2.add(peak6)
    peaks2.add(peak8)
    peaks1, peaks2 = classify_peaks_by_overlap(peaks1, peaks2)
    peaks_merged = merge_common_peaks(peaks1, peaks2)
    assert peaks_merged.size == 2
    assert peaks_merged.fetch('chr1')[0].start == 1
    assert peaks_merged.fetch('chr1')[0].end == 200
    assert peaks_merged.fetch('chr1')[0].summit == 143
    assert peaks_merged.fetch('chr2')[0].start == 1
    assert peaks_merged.fetch('chr2')[0].end == 120
    assert peaks_merged.fetch('chr2')[0].summit == 65


def test_generate_random_peaks():
    peaks_ref = GenomicRegions(name='test1')
    peak1 = ManormPeak(chrom='chr1', start=1, end=100, summit=50)
    peak2 = ManormPeak(chrom='chr1', start=50, end=150, summit=120)
    peak3 = ManormPeak(chrom='chr1', start=149, end=200, summit=166)
    peak4 = ManormPeak(chrom='chr1', start=250, end=500, summit=350)
    peak5 = ManormPeak(chrom='chr2', start=300, end=400, summit=333)
    peak6 = ManormPeak(chrom='chr2', start=1, end=100, summit=50)
    peak7 = ManormPeak(chrom='chr2', start=50, end=120, summit=80)
    peak8 = ManormPeak(chrom='chr22', start=1, end=100, summit=50)
    peaks_ref.add(peak1)
    peaks_ref.add(peak2)
    peaks_ref.add(peak3)
    peaks_ref.add(peak4)
    peaks_ref.add(peak5)
    peaks_ref.add(peak6)
    peaks_ref.add(peak7)
    peaks_ref.add(peak8)
    peaks_rand = generate_random_regions(peaks_ref)
    assert peaks_rand.size == 8
    assert len(peaks_rand.fetch('chr1')) == 4
    assert len(peaks_rand.fetch('chr2')) == 3
    assert len(peaks_rand.fetch('chr22')) == 1


def test_peaks_count_n_unique_common():
    peaks = GenomicRegions(name='test')
    peak1 = ManormPeak(chrom='chr1', start=1, end=100)
    peak1.iscommon = False
    peak2 = ManormPeak(chrom='chr1', start=10, end=50)
    peak2.iscommon = True
    peaks.add(peak1)
    peaks.add(peak2)
    assert count_unique_peaks(peaks) == 1
    assert count_common_peaks(peaks) == 1
    peak3 = ManormPeak(chrom='chr11', start=100, end=200)
    peak3.iscommon = False
    peaks.add(peak3)
    assert count_unique_peaks(peaks) == 2
    assert count_common_peaks(peaks) == 1
    peak4 = ManormPeak(chrom='chr12', start=100, end=200)
    peak4.iscommon = True
    peaks.add(peak4)
    assert count_unique_peaks(peaks) == 2
    assert count_common_peaks(peaks) == 2
