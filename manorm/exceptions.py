# -*- coding: utf-8 -*-

"""
manorm.exceptions
~~~~~~~~~~~~~~~~~

This module contains the exceptions of MAnorm.
"""


class PeakFormatError(Exception):
    """Invalid peak format error."""

    def __init__(self, format, line):
        self.format = format
        self.line = line.rstrip()

    def __str__(self):
        return "Invalid {} format: {!r}".format(self.format, self.line)


class UnknownFormatError(Exception):
    """Unknown format error."""

    def __init__(self, format):
        self.format = format

    def __str__(self):
        return "Unknown format: {!r}".format(self.format)
