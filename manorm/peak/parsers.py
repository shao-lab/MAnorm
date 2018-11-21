"""
manorm.peak.parsers
~~~~~~~~~~~~~~~~~~~

Peak parsers to read/write peak files.
"""
from __future__ import absolute_import, division

import gzip
import logging
import os

from manorm import compat
from manorm.exceptions import PeakFormatError, UnknownFormatError
from manorm.peak import Peak

logger = logging.getLogger(__name__)


def get_peak_parser(path, format):
    try:
        return peak_parsers[format](path)
    except KeyError:
        raise UnknownFormatError(format=format)


class PeakParser(object):
    """Peak parser base class."""

    def __init__(self, path):
        self.path = path
        self.compressed = True
        self.field_count = 0
        with gzip.open(path, 'rb') as fin:
            try:
                fin.read(1)
                logger.debug("Input file {} is gzipped".format(os.path.basename(self.path)))
            except IOError:
                self.compressed = False
        if self.compressed:
            self.handle = gzip.open(self.path, 'rb')
        else:
            self.handle = compat.open(self.path, 'rb')

    def close(self):
        self.handle.close()

    @staticmethod
    def _is_header(line):
        """Check whether the given line is a header line."""

        if line == '':  # skip empty lines
            return True
        if line.startswith('#') or line.startswith('track') or line.startswith('browser'):  # BED custom track header
            return True
        if line.strip().split('\t')[0] == 'chr':  # MACS-xls header
            return True
        return False

    def _parse_line(self, line):
        """Abstract method to parse single line into a peak."""
        raise NotImplementedError

    def parse(self, rows=0):
        """Parse lines to get peaks from the input peak file.

        Args:
            rows (int): Only parse top {rows} lines. Set it to 0 (default) to parse all lines.

        Returns:
            peaks (list): [peak1, peak2, ...]

        """

        header = True
        peaks = []
        self.handle.seek(0)  # parse from the head of the file
        line_num = 0
        for line in self.handle:
            line_num += 1

            if header and self._is_header(line):
                logger.debug("Detected header line: [{}]".format(line.rstrip()))
                continue
            if header:
                logger.debug("Detected first data entry: [{}]".format(line.rstrip()))
                self.field_count = len(line.strip().split('\t'))
                logger.debug("Number of columns: {}".format(self.field_count))
                header = False
            peaks.append(self._parse_line(line))
            if 0 < rows <= line_num:  # only parse top {rows} lines
                break
        return peaks

    def check(self, rows=0):
        """Check whether the input peak file is matched with corresponding format."""
        try:
            self.parse(rows)
        except PeakFormatError:
            return False
        self.close()
        return True


class BEDParser(PeakParser):
    """Peak parser for BED format."""

    def _parse_line(self, line):
        fields = line.strip().split('\t')  # BED required fields to be tab-delimited
        field_count = len(fields)
        try:
            assert field_count >= 3  # BED requires at least 3 columns for chrom, start, end
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = '.'
            score = 0
            strand = '.'

            if field_count > 3:
                name = fields[3]
                if field_count > 4:
                    score = float(fields[4])
                    if field_count > 5:
                        strand = fields[5]

            peak = Peak(chrom=chrom, start=start, end=end, name=name, score=score, strand=strand)
        except:
            self.close()
            raise PeakFormatError(line)
        return peak


class BEDSummitParser(PeakParser):
    """Peak parser for BED-summit format."""

    def _parse_line(self, line):
        fields = line.strip().split('\t')
        field_count = len(fields)
        try:
            assert field_count >= 4  # BED-summit requires at least 4 columns for chrom, start, end, summit
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            summit = int(fields[3]) + start
            name = '.'
            score = 0
            strand = '.'
            if field_count > 4:
                name = fields[4]
                if field_count > 5:
                    score = float(fields[5])
                    if field_count > 6:
                        strand = fields[6]
            peak = Peak(chrom=chrom, start=start, end=end, summit=summit, name=name, score=score, strand=strand)
        except:
            self.close()
            raise PeakFormatError(line)
        return peak


class MACSParser(PeakParser):
    """Peak parser for MACS-xls format."""

    def _parse_line(self, line):
        fields = line.strip().split('\t')
        field_count = len(fields)
        try:
            assert field_count >= 5  # MACS-xls requires at least 5 columns for chrom, start, end, length, summit
            chrom = fields[0]
            start = int(fields[1]) - 1  # 1-based coordinates
            end = int(fields[2])
            summit = int(fields[4]) + start
            name = '.'
            score = 0
            strand = '.'
            if field_count > 6:
                score = float(fields[6])  # score represents the -10*log10(p-value) of macs peaks
            peak = Peak(chrom=chrom, start=start, end=end, summit=summit, name=name, score=score, strand=strand)
        except:
            self.close()
            raise PeakFormatError(line)
        return peak


class MACS2Parser(PeakParser):
    pass


class NarrowPeakParser(PeakParser):
    pass


class BroadPeakParser(PeakParser):
    pass


peak_parsers = {'bed': BEDParser,
                'bed-summit': BEDSummitParser,
                'macs': MACSParser,
                'macs2': MACS2Parser,
                'narrowpeak': NarrowPeakParser,
                'broadpeak': BroadPeakParser}
