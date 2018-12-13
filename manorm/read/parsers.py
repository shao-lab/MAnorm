# -*- coding: utf-8 -*-

"""
manorm.read.parsers
~~~~~~~~~~~~~~~~~~~

Read parsers to read peak files.
"""


from __future__ import absolute_import, division

import pysam
import gzip
import logging
import os

from manorm.compat import open
from manorm.exceptions import InvalidFormatError

logger = logging.getLogger(__name__)


class BEDParser(object):
    """Class for read BED files."""

    def __init__(self, path):
        self.path = path
        self.format = 'BED'
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
        elif line.startswith('#') or line.startswith('track') or line.startswith('browser'):  # BED header
            return 2
        return 0

    @staticmethod
    def _parse_line(line, shift=100, *args, **kwargs):
        fields = line.rstrip().split('\t')  # BED required fields to be tab-delimited
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[5]
        if strand == '+':
            pos = start + shift
        elif strand == '-':
            pos = end - shift
        else:
            raise ValueError
        return chrom, pos

    def parse(self, *args, **kwargs):
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
                yield self._parse_line(line=line, *args, **kwargs)
            except (IndexError, ValueError):
                self.close()
                raise InvalidFormatError(format=self.format, line=line, line_num=line_num)


class BEDPEParser(BEDParser):
    """Class for read BEDPE files."""

    def __init__(self, path):
        super(BEDPEParser, self).__init__(path)
        self.format = 'BEDPE'

    @staticmethod
    def _parse_line(line, *args, **kwargs):
        fields = line.rstrip().split('\t')  # BEDPE required fields to be tab-delimited
        chrom1 = fields[0]
        start1 = int(fields[1])
        end1 = int(fields[2])
        chrom2 = fields[3]
        start2 = int(fields[4])
        end2 = int(fields[5])
        start = min(start1, start2)
        end = max(end1, end2)
        if chrom1 == chrom2:
            return chrom1, (start+end) // 2
        else:
            return None, None


class SAMParser(object):
    pass

class BAMParser(object):
    pass