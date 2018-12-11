# -*- coding: utf-8 -*-

"""
manorm.exceptions
~~~~~~~~~~~~~~~~~

This module contains the exceptions of MAnorm.
"""


class PeakFormatError(Exception):
    """Invalid peak format error."""

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


class UnknownFormatError(Exception):
    """Unknown format error."""

    def __init__(self, format):
        self.format = format

    def __str__(self):
        return "Unknown format: {!r}".format(self.format)
