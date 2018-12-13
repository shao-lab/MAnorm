# -*- coding: utf-8 -*-

"""
manorm.exceptions
~~~~~~~~~~~~~~~~~

This module contains the exceptions of MAnorm.
"""


class InvalidFormatError(Exception):
    """Invalid file format error."""

    def __init__(self, format, line, line_num=None):
        self.format = format
        self.line = line.rstrip()
        self.line_num = line_num

    def __str__(self):
        if self.line_num is None:
            msg = "Invalid {} format: {!r}".format(self.format, self.line)
        else:
            msg = "Invalid {} format at line {}: {!r}".format(self.format, self.line_num, self.line)
        return msg


class UnsupportedFormatError(Exception):
    """Unsupported format error."""

    def __init__(self, format):
        self.format = format

    def __str__(self):
        return "Unsupported format: {!r}".format(self.format)


class UnmatchedBedFormatError(Exception):
    """Unmatched format error."""

    def __str__(self):
        return "Unmatched format: please use BED format for single-end reads and BEDPE for paired-end reads"
