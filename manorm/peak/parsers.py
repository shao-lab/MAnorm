# -*- coding: utf-8 -*-

"""
manorm.peak.parsers
~~~~~~~~~~~~~~~~~~~

Parsers to parse peak files.
"""

from __future__ import absolute_import

import gzip
import logging
import os

from manorm.compat import open
from manorm.exceptions import InvalidFormatError

logger = logging.getLogger(__name__)


class PeakParser(object):
    """Base class for peak file parsers."""

    def __init__(self, path):
        self.path = path
        self.format = None
        self.compressed = True
        with gzip.open(path, 'rb') as fin:
            try:
                fin.read(1)
                logger.debug("Input file is gzipped: {}".format(os.path.abspath(self.path)))
            except IOError:
                self.compressed = False
        if self.compressed:
            self.handle = gzip.open(self.path, 'rb')
        else:
            self.handle = open(self.path)

    def close(self):
        """Close the peak file."""
        self.handle.close()

    @staticmethod
    def _is_header(line):
        """Check if a line is header."""
        line = line.rstrip()
        if line == '':  # skip empty lines
            return 1
        elif (line.startswith('#') or line.startswith('track') or line.startswith('browser')  # BED header
              or line.split('\t')[0] == 'chr'):  # MACS/MACS2 xls header
            return 2
        return 0

    @staticmethod
    def _parse_line(line):
        """Abstract method to parse a line."""
        raise NotImplementedError

    def parse(self):
        """Parse lines to get peaks from the input peak file."""
        self.handle.seek(0)  # parse from the beginning of the file
        line_num = 0
        expect_header = True
        for line in self.handle:
            line_num += 1
            if expect_header:
                status = self._is_header(line)
                if status == 1:
                    logger.debug("Detected empty line at line {}".format(line_num))
                    continue
                elif status == 2:
                    logger.debug("Detected header line at line {}: {}".format(line_num, line.rstrip()))
                    continue
                else:
                    expect_header = False
            try:
                yield self._parse_line(line=line)
            except (IndexError, ValueError):
                self.close()
                raise InvalidFormatError(format=self.format, line=line, line_num=line_num)


class BEDParser(PeakParser):
    """Peak parser for BED format."""

    def __init__(self, path):
        super(BEDParser, self).__init__(path)
        self.format = 'BED'

    @staticmethod
    def _parse_line(line):
        fields = line.rstrip().split('\t')  # BED required fields to be tab-delimited
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        summit = None
        return chrom, start, end, summit


class BEDSummitParser(PeakParser):
    """Peak parser for BED-summit format."""

    def __init__(self, path):
        super(BEDSummitParser, self).__init__(path)
        self.format = 'BED-summit'

    @staticmethod
    def _parse_line(line):
        fields = line.rstrip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        summit = int(fields[3]) + start
        return chrom, start, end, summit


class MACSParser(PeakParser):
    """Peak parser for MACS-xls format."""

    def __init__(self, path):
        super(MACSParser, self).__init__(path)
        self.format = 'MACS'

    @staticmethod
    def _parse_line(line):
        fields = line.rstrip().split('\t')
        chrom = fields[0]
        start = int(fields[1]) - 1  # 1-based coordinates
        end = int(fields[2])
        summit = int(fields[4]) + start  # relative peak summit position (from MACS manual)
        return chrom, start, end, summit


class MACS2Parser(PeakParser):
    """Peak parser for MACS2-xls format."""

    def __init__(self, path):
        super(MACS2Parser, self).__init__(path)
        self.format = 'MACS2'

    @staticmethod
    def _parse_line(line):
        fields = line.rstrip().split('\t')
        chrom = fields[0]
        start = int(fields[1]) - 1  # 1-based coordinates
        end = int(fields[2])
        summit = int(fields[4])  # absolute peak summit position (from MACS2 manual)
        return chrom, start, end, summit


class NarrowPeakParser(PeakParser):
    """Peak parser for NarrowPeak format."""

    def __init__(self, path):
        super(NarrowPeakParser, self).__init__(path)
        self.format = 'NarrowPeak'

    @staticmethod
    def _parse_line(line):
        fields = line.rstrip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        summit = int(fields[9]) + start  # https://genome.ucsc.edu/FAQ/FAQformat.html#format12
        if summit == start - 1:
            summit = None
        return chrom, start, end, summit
