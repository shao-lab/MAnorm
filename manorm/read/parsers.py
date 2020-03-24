"""
manorm.read.parsers
-------------------

Parsers to parse the sequencing reads file.
"""

import gzip
import logging

import pysam

from manorm.exceptions import FileFormatError

logger = logging.getLogger(__name__)


class BedReadParser:
    """Read parser for the BED format."""

    def __init__(self, path):
        self.path = path
        self.format = 'BED'

    @property
    def is_gzipped(self):
        with gzip.open(self.path, 'rb') as fin:
            try:
                fin.read(1)
                is_gzipped = True
            except OSError:
                is_gzipped = False
        return is_gzipped

    @staticmethod
    def _is_header(line):
        """Check if a line is header."""
        line = line.strip()
        if line.startswith('#') or line.startswith('track') or line.startswith(
                'browser'):  # BED header
            return True
        else:
            return False

    @staticmethod
    def _parse_line(line, shift=100, *args, **kwargs):
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[5]
        if strand == '+':
            pos = start + shift
        elif strand == '-':
            pos = end - shift
        else:
            raise ValueError(f"invalid strand: {strand!r}")
        return chrom, pos

    def parse(self, *args, **kwargs):
        """Parse lines to get reads from the input read file."""
        if self.is_gzipped:
            fin = gzip.open(self.path, 'rb')
        else:
            fin = open(self.path, 'r')
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
                chrom, pos = self._parse_line(line, *args, **kwargs)
                if chrom is not None:
                    yield chrom, pos
            except (IndexError, ValueError, TypeError):
                fin.close()
                raise FileFormatError(format=self.format, line_num=line_num,
                                      line=line)
        fin.close()


class BedPeReadParser(BedReadParser):
    """Read parser for the BEDPE format."""

    def __init__(self, path):
        super().__init__(path)
        self.format = 'BEDPE'

    @staticmethod
    def _parse_line(line, *args, **kwargs):
        fields = line.strip().split('\t')
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


class SamReadParser:
    """Read parser for the SAM format.
    Ref: http://samtools.sourceforge.net/SAM1.pdf
    """

    def __init__(self, path):
        self.path = path
        self.format = 'SAM'
        self.handle = pysam.AlignmentFile(self.path, 'r')

    def parse(self, paired=False, shift=100):
        """Parse lines to get reads from the input read file."""
        for read in self.handle:
            if read.is_unmapped or read.is_qcfail or read.is_secondary \
                    or read.is_supplementary:
                continue
            if paired:
                if not read.is_paired:
                    logger.debug(
                        f"Skipped single-end read: {read.to_string()!r}")
                    continue
                if read.is_read1 and read.is_proper_pair \
                        and not read.mate_is_unmapped:
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
                    if read.template_length == 0:
                        logger.debug(
                            f"Detected read with TLEN=0: {read.to_string()!r}")
                    if chrom1 == chrom2:
                        yield chrom1, (start + end) // 2
                    else:
                        continue
                else:
                    continue
            else:
                if read.is_paired:
                    logger.debug(
                        f"Skipped paired-end read: {read.to_string()!r}")
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
        self.handle.close()


class BamReadParser(SamReadParser):
    """Read parser for BAM format."""

    def __init__(self, path):
        self.path = path
        self.format = 'BAM'
        self.handle = pysam.AlignmentFile(self.path, 'rb')


def get_read_parser(format):
    """Get proper read parser for the given format.

    Parameters
    ----------
    format : str
        File format (case-insensitive).

    Returns
    -------
        Corresponding read parser.
    """
    format = format.lower()
    if format == 'bed':
        return BedReadParser
    elif format == 'bedpe':
        return BedPeReadParser
    elif format == 'sam':
        return SamReadParser
    elif format == 'bam':
        return BamReadParser
    else:
        raise ValueError(f"unknown read file format: {format!r}")
