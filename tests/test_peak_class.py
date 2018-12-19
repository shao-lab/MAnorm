from manorm.peak import Peak, Peaks


def test_peak_init():
    peak = Peak(chrom='chr1', start=1, end=100)
    assert peak.summit == 50


def test_peaks_init():
    peaks = Peaks(name='test')
    assert peaks.name == 'test'
    assert peaks.size == 0
    assert peaks.chroms == []


def test_peaks_add():
    peaks = Peaks(name='test')
    peak = Peak(chrom='chr1', start=1, end=100)
    peaks.add(peak)
    assert peaks.size == 1
    assert peaks.chroms == ['chr1']
    assert peaks.data['chr1'] == [peak]


def test_peaks_fetch():
    peaks = Peaks(name='test')
    peak = Peak(chrom='chr1', start=1, end=100)
    peaks.add(peak)
    assert peaks.size == 1
    assert peaks.chroms == ['chr1']
    assert peaks.fetch(chrom='chr1') == [peak]
    assert peaks.fetch(chrom='chr11') == []


def test_peaks_sort():
    peaks = Peaks(name='test')
    peak1 = Peak(chrom='chr1', start=1, end=100)
    peak2 = Peak(chrom='chr1', start=10, end=50)
    peaks.add(peak=peak1)
    peaks.add(peak=peak2)
    assert peaks.size == 2
    assert peaks.chroms == ['chr1']
    peaks.sort(by='start', ascending=True)
    assert peaks.fetch(chrom='chr1')[0] == peak1
    peaks.sort(by='start', ascending=False)
    assert peaks.fetch(chrom='chr1')[0] == peak2
    peaks.sort(by='end', ascending=True)
    assert peaks.fetch(chrom='chr1')[0] == peak2
    peaks.sort(by='end', ascending=False)
    assert peaks.fetch(chrom='chr1')[0] == peak1
    peaks.sort(by='summit', ascending=True)
    assert peaks.fetch(chrom='chr1')[0] == peak2
    peaks.sort(by='summit', ascending=False)
    assert peaks.fetch(chrom='chr1')[0] == peak1


def test_peaks_n_unique():
    peaks = Peaks(name='test')
    peak1 = Peak(chrom='chr1', start=1, end=100)
    peak1.type = 'unique'
    peak2 = Peak(chrom='chr1', start=10, end=50)
    peak2.type = 'common'
    peaks.add(peak=peak1)
    peaks.add(peak=peak2)
    assert peaks.n_unique == 1
    assert peaks.n_common == 1
    peak3 = Peak(chrom='chr11', start=100, end=200)
    peak3.type = 'unique'
    peaks.add(peak=peak3)
    assert peaks.n_unique == 2
    assert peaks.n_common == 1
    peak4 = Peak(chrom='chr12', start=100, end=200)
    peak4.type = 'common'
    peaks.add(peak=peak4)
    assert peaks.n_unique == 2
    assert peaks.n_common == 2
