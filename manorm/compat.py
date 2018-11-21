"""Code compatibility for Python 2 and 3."""

from __future__ import absolute_import

import functools
import sys
import itertools

PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3

if PY3:
    long = int

    text_type = str
    string_types = (str,)
    integer_types = (int,)

    range = range
    map = map
    zip = zip
    filter = filter
    reduce = functools.reduce

else:
    # Python 2
    text_type = unicode
    string_types = (str, unicode)
    integer_types = (int, long)

    range = xrange
    map = itertools.imap
    zip = itertools.izip
    filter = itertools.ifilter
    reduce = reduce
