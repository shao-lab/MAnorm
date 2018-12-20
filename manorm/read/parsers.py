# -*- coding: utf-8 -*-

"""
manorm.read.parsers
~~~~~~~~~~~~~~~~~~~

Parsers to parse read files.
"""

from __future__ import absolute_import, division

import gzip
import logging
import os

import pysam

from manorm.compat import open
from manorm.exceptions import InvalidFormatError

logger = logging.getLogger(__name__)


class BEDParser(object):
    """Read parser for BED format."""

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
        """Close the read file."""
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
        """Parse lines to get reads from the input read file."""
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
                chrom, pos = self._parse_line(line=line, *args, **kwargs)
                if chrom is not None:
                    yield chrom, pos
            except (IndexError, ValueError):
                self.close()
                raise InvalidFormatError(format=self.format, line=line, line_num=line_num)


class BEDPEParser(BEDParser):
    """Read parser for BEDPE format."""

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
            return chrom1, (start + end) // 2
        else:
            return None, None


class SAMParser(object):
    """Read parser for SAM format.
    Ref: http://samtools.sourceforge.net/SAM1.pdf
    """

    def __init__(self, path):
        self.path = path
        self.format = 'SAM'
        self.handle = pysam.AlignmentFile(self.path, 'r')

    def close(self):
        """Close the read file."""
        self.handle.close()

    def parse(self, paired=False, shift=100):
        """Parse lines to get reads from the input read file."""
        for read in self.handle:
            if read.is_unmapped or read.is_qcfail or read.is_secondary or read.is_supplementary:
                continue
            if paired:
                if not read.is_paired:
                    logger.warning(
                        "Detected single-end read in paired-end mode: {!r}, skipped.".format(read.to_string()))
                    continue
                if read.is_read1 and read.is_proper_pair and not read.mate_is_unmapped:
                    chrom1 = read.reference_name
                    start1 = read.reference_start
                    end1 = read.reference_end
                    chrom2 = read.next_reference_name
                    if read.is_reverse:
                        start = end1 + read.template_length
                        end = end1
                    else:
                        start = start1
                        end = start1 + read.template_length
                    if chrom1 == chrom2:
                        yield chrom1, (start + end) // 2
                    else:
                        continue
                else:
                    continue
            else:
                if read.is_paired:
                    logger.warning(
                        "Detected paired-end read in single-end mode: {!r}, skipped.".format(read.to_string()))
                    continue
                if read.is_unmapped:
                    continue
                else:
                    chrom = read.reference_name
                    start = read.reference_start
                    end = read.reference_end
                    if read.is_reverse:
                        pos = end - shift
                    else:
                        pos = start + shift
                    yield chrom, pos


class BAMParser(SAMParser):
    """Read parser for BAM format."""

    def __init__(self, path):
        self.path = path
        self.format = 'BAM'
        self.handle = pysam.AlignmentFile(self.path, 'rb')
