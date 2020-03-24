"""
manorm.region
-------------

Module for genomic regions and MAnorm peaks.
"""

import logging
import os

from manorm.region.parsers import get_region_parser
from manorm.stats import xy_to_ma, ma_to_xy, manorm_p

logger = logging.getLogger(__name__)

REGION_FORMATS = ['bed', 'bed3-summit', 'macs', 'macs2', 'narrowpeak',
                  'broadpeak']


class GenomicRegion:
    """Class for a genomic region.

    Parameters
    ----------
    chrom : str
        The chromosome name of the region.
    start : int
        The start coordinate of the region.
    end : int
        The end coordinate of the region.
    summit : int, optional
        The summit coordinate of the region. If not specified, the middle point
        will be taken as the summit.

    Attributes
    ----------
    chrom : str
        The chromosome name of the region.
    start : int
        The start coordinate of the region.
    end : int
        The end coordinate of the region.
    summit : int
        The summit coordinate of the region.

    Notes
    -----
    The coordinates are 0-based, which means the region range is [start, end)
    and `start` <= `summit` < `end` is strictly inspected.
    """

    def __init__(self, chrom, start, end, summit=None):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        if self.start >= self.end:
            raise ValueError(
                f"expect start < end, got: start={start} end={end}")
        if summit is not None:
            self.summit = int(summit)
        else:
            self.summit = (self.start + self.end) // 2
        if not self.start <= self.summit < self.end:
            raise ValueError(
                f"expect start <= summit < end, got start={start} "
                f"summit={summit} end={end}")

    def __repr__(self):
        return f"GenomicRegion({self.chrom}:{self.start}-{self.end})"


class ManormPeak(GenomicRegion):
    """Class for a MAnorm peak.

    Parameters
    ----------
    chrom : str
        The chromosome name of the peak.
    start : int
        The start coordinate of the peak.
    end : int
        The end coordinate of the peak.
    summit : int, optional
        The summit coordinate of the peak. If not specified, the middle point
        will be taken as the summit.

    Attributes
    ----------
    chrom : str
        The chromosome name of the peak.
    start : int
        The start coordinate of the peak.
    end : int
        The end coordinate of the peak.
    summit : int
        The summit coordinate of the peak.
    counted : bool
        Whether the reads are counted or not.
    read_count1 : int or None
        The number of reads falling into the peak region in sample 1.
    read_count2 : int or None
        The number of reads falling into the peak region in sample 2.
    read_density1 : float or None
        The read density of sample 1 in the unit of reads per kilobase.
    read_density2 : float or None
        The read density of sample 2 in the unit of reads per kilobase.
    m_raw: float or None
        The raw(unnormalized) M value.
    a_raw: float or None
        The raw(unnormalized) A value.
    normalized : bool
        Whether the peak is normalized.
    read_density1_normed : float or None
        Normalized read density of sample 1.
    read_density2_normed : float or None
        Normalized read density of sample 2.
    m_normed : float or None
        Normalized M value.
    a_normed : float or None
        Normalized A value.
    p_value : float or None
        The P value of MAnorm.
    iscommon : bool
        Indicator to show whether the peak is a common peak.
    summit_dis : int or None
        The minimal distance between the summits of two common peaks, only
        valid for merged common peaks.

    Notes
    -----
    After counting the reads located in the peak region, fields `read_count1`,
    `read_count2`, `read_density1`, `read_density2`, `m_raw`, `a_raw` are set
    to proper values.
    After the normalization step, fields `read_density1_normed`,
    `read_density2_normed`, `m_normed`, `a_normed`, `p_value` are set to
    proper values.
    """

    def __init__(self, chrom, start, end, summit=None):
        super().__init__(chrom, start, end, summit)
        self.counted = False
        self.read_count1 = None
        self.read_count2 = None
        self.read_density1 = None
        self.read_density2 = None
        self.m_raw = None
        self.a_raw = None
        self.normalized = False
        self.read_density1_normed = None
        self.read_density2_normed = None
        self.m_normed = None
        self.a_normed = None
        self.p_value = None
        self.iscommon = False
        self.summit_dis = None

    def count_reads(self, reads1, reads2, window=2000):
        """Count reads, calculate the read densities and raw (M, A) values.

        Parameters
        ----------
        reads1 : `manorm.read.Reads`
            MAnorm reads object of sample 1.
        reads2 : `manorm.read.Reads`
            MAnorm reads object of sample 2.
        window : int, optional
            The window size to count reads, default=2000.
        """
        if window <= 0:
            raise ValueError(f"expect window size > 0, got {window}")
        extend = window // 2
        self.read_count1 = reads1.count(self.chrom, self.summit - extend,
                                        self.summit + extend) + 1
        self.read_count2 = reads2.count(self.chrom, self.summit - extend,
                                        self.summit + extend) + 1
        self.read_density1 = self.read_count1 * 1000 / (extend * 2)
        self.read_density2 = self.read_count2 * 1000 / (extend * 2)
        self.m_raw, self.a_raw = xy_to_ma(self.read_density1,
                                          self.read_density2)
        self.counted = True

    def normalize(self, slope, intercept):
        """Normalize the M value and A value to remove the global dependence of
        M on A given the coefficients of fitted robust linear model.

        M-A model: M = slope * A + intercept.

        Parameters
        ----------
        slope : float
            The slope of fitted M-A linear model.
        intercept : float
            The intercept of fitted M-A linear model.
        """
        self.m_normed = self.m_raw - (slope * self.a_raw + intercept)
        self.a_normed = self.a_raw
        self.read_density1_normed, self.read_density2_normed = ma_to_xy(
            self.m_normed, self.a_normed)
        self.p_value = manorm_p(self.read_density1_normed,
                                self.read_density2_normed)
        self.normalized = True


class GenomicRegions:
    """Class for a collection of genomic regions.

    Parameters
    ----------
    name : str, optional
        The name of the genomic regions.

    Attributes
    ----------
    name : str or None
        The name of the genomic regions.
    """

    def __init__(self, name=None):
        self.name = name
        self._data = {}

    @property
    def chroms(self):
        """Returns sorted chromosome names of the genomic regions.

        Returns
        -------
        list of str
            Chromosome names (sorted) of the genomic regions.
        """
        return sorted(self._data.keys())

    @property
    def size(self):
        """Returns the number of genomic regions."""
        return sum(len(value) for value in self._data.values())

    def add(self, region):
        """Add a genomic region into the collection.

        Parameters
        ----------
        region : GenomicRegion
            The genomic region to be added into the collection.
        """
        if not isinstance(region, GenomicRegion):
            raise ValueError("requires a `GenomicRegion` object to be added")
        else:
            self._data.setdefault(region.chrom, [])
            self._data[region.chrom].append(region)

    def sort(self, by='start', ascending=True):
        """Sort genomic regions.

        Parameters
        ----------
        by : str, optional
            Which atrribute is used to sort by, default='start'.
        ascending : bool, optional
            Sort ascendingly or not, default=True.
        """
        for chrom in self.chroms:
            self._data[chrom].sort(key=lambda x: getattr(x, by),
                                   reverse=not ascending)

    def fetch(self, chrom):
        """Fetch genomic regions on specified chromosome.

        Parameters
        ----------
        chrom : str
             The chromosome name to fetch regions from.

        Returns
        -------
        list
            A list of genomic regions on the specified chromosome.
        """
        if chrom in self._data:
            return self._data[chrom]
        else:
            return []


def load_genomic_regions(path, format='bed', name=None):
    """Read genomic regions from the specified path.

    Parameters
    ----------
    path : str
        Path to load the genomic regions.
    format : str, optional
        File format, default='bed'.
    name : str, optional
        Name of the genomic regions. If not specified, the basename of the
        file will be used.

    Returns
    -------
    regions : GenomicRegions
        Loaded genomic regions.
    """
    logger.info(f"Loading genomic regions from {path} [{format}]")
    if name is None:
        name = os.path.splitext(os.path.basename(path))[0]
    parser = get_region_parser(format)()
    regions = GenomicRegions(name)
    for chrom, start, end, summit in parser.parse(path):
        regions.add(GenomicRegion(chrom, start, end, summit))
    regions.sort()
    logger.info(f"Loaded {regions.size} genomic regions")
    return regions


def load_manorm_peaks(path, format='bed', name=None):
    """Read peaks from the specified path.

    Parameters
    ----------
    path : str
        Path to load the peaks.
    format : str, optional
        File format, default='bed'.
    name : str, optional
        Name of the peaks. If not specified, the basename of the file will be
        used.

    Returns
    -------
    peaks : GenomicRegions
        Loaded peaks.
    """
    logger.info(f"Loading peaks from {path} [{format}]")
    if name is None:
        name = os.path.splitext(os.path.basename(path))[0]
    parser = get_region_parser(format)()
    peaks = GenomicRegions(name)
    for chrom, start, end, summit in parser.parse(path):
        peaks.add(ManormPeak(chrom, start, end, summit))
    peaks.sort()
    logger.info(f"Loaded {peaks.size} peaks")
    return peaks
