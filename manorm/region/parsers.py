"""
manorm.region.parsers
---------------------

Parsers to read the genomic region files.
"""

import logging

from manorm.exceptions import FileFormatError

logger = logging.getLogger(__name__)


def is_track_header(line):
    """Returns if the line is a header line used in genome tracks/browers."""
    line = line.strip()
    if line.startswith('#') or line.startswith('track') or line.startswith(
            'browser'):
        return True
    else:
        return False


def is_comment_header(line):
    """Returns if the line is a comment header line."""
    line = line.strip()
    if line.startswith('#'):
        return True
    else:
        return False


def is_macs_header(line):
    """Returns if the line is a header line used in MACS/MACS2 xls."""
    line = line.strip()
    if line.startswith('#') or line.split('\t')[0] == 'chr':
        return True
    else:
        return False


class RegionParser:
    """Base class for region file parsers."""

    def __init__(self, format):
        self.format = format

    @staticmethod
    def _is_header(line):
        """Abstract method to check if a line is a header line."""
        raise NotImplementedError

    @staticmethod
    def _parse_line(line):
        """Abstract method to parse a line."""
        raise NotImplementedError

    def parse(self, path):
        """Read genomic regions from the given file."""
        with open(path, 'r') as fin:
            line_num = 0
            expect_header = True
            for line in fin:
                line_num += 1
                line = line.strip()
                if not line:  # skip empty lines
                    continue
                if expect_header:
                    if self._is_header(line):
                        logger.debug(
                            f"Detected header at line {line_num}: {line!r}")
                        continue
                    else:
                        expect_header = False
                try:
                    yield self._parse_line(line)
                except (IndexError, ValueError, TypeError):
                    raise FileFormatError(format=self.format,
                                          line_num=line_num, line=line)


class BedRegionParser(RegionParser):
    """Region parser for the BED format."""

    def __init__(self):
        super().__init__('BED')

    @staticmethod
    def _is_header(line):
        return is_track_header(line)

    @staticmethod
    def _parse_line(line):
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        summit = None
        return chrom, start, end, summit


class Bed3SummitRegionParser(RegionParser):
    """Region parser for the BED3-summit format."""

    def __init__(self):
        super().__init__(format='BED3-summit')

    @staticmethod
    def _is_header(line):
        return is_comment_header(line)

    @staticmethod
    def _parse_line(line):
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        summit = int(fields[3])
        return chrom, start, end, summit


class MacsRegionParser(RegionParser):
    """Region parser for the MACS-xls format."""

    def __init__(self):
        super().__init__(format='MACS-xls')

    @staticmethod
    def _is_header(line):
        return is_macs_header(line)

    @staticmethod
    def _parse_line(line):
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1]) - 1  # coordinates are 1-based in MACS xls
        end = int(fields[2])
        summit = int(fields[4]) + start  # relative summit pos for MACS1
        return chrom, start, end, summit


class Macs2RegionParser(RegionParser):
    """Region parser for the MACS2-xls format."""

    def __init__(self):
        super().__init__(format='MACS2-xls')

    @staticmethod
    def _is_header(line):
        return is_macs_header(line)

    @staticmethod
    def _parse_line(line):
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1]) - 1  # coordinates are 1-based in MACS2 xls
        end = int(fields[2])
        summit = int(fields[4]) - 1  # absolute summit pos for MACS2
        return chrom, start, end, summit


class NarrowPeakRegionParser(RegionParser):
    """Region parser for the NarrowPeak format."""

    def __init__(self):
        super().__init__(format='NarrowPeak')

    @staticmethod
    def _is_header(line):
        return is_track_header(line)

    @staticmethod
    def _parse_line(line):
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        # https://genome.ucsc.edu/FAQ/FAQformat.html#format12
        summit = int(fields[9])
        if summit == -1:
            summit = None
        else:
            summit = start + summit
        return chrom, start, end, summit


class BroadPeakRegionParser(RegionParser):
    """Region parser for the BroadPeak format."""

    def __init__(self):
        super().__init__(format='BroadPeak')

    @staticmethod
    def _is_header(line):
        return is_track_header(line)

    @staticmethod
    def _parse_line(line):
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        # https://genome.ucsc.edu/FAQ/FAQformat.html#format13
        summit = None
        return chrom, start, end, summit


def get_region_parser(format):
    """Get proper region parser for the given format.

    Parameters
    ----------
    format : str
        File format (case-insensitive).

    Returns
    -------
        Corresponding region parser.
    """
    format = format.lower()
    if format == 'bed':
        return BedRegionParser
    elif format == 'bed3-summit':
        return Bed3SummitRegionParser
    elif format == 'macs':
        return MacsRegionParser
    elif format == 'macs2':
        return Macs2RegionParser
    elif format == 'narrowpeak':
        return NarrowPeakRegionParser
    elif format == 'broadpeak':
        return BroadPeakRegionParser
    else:
        raise ValueError(f"unknown region file format: {format!r}")
