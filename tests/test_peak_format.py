import os

import pytest

from manorm.exceptions import InvalidFormatError, UnsupportedFormatError
from manorm.peak import load_peaks


def test_unsupported_format(data_dir):
    with pytest.raises(UnsupportedFormatError):
        peaks = load_peaks(os.path.join(data_dir, 'test_peaks.bed'), format='unknown_format')


def test_invalid_format(data_dir):
    with pytest.raises(InvalidFormatError):
        peaks = load_peaks(os.path.join(data_dir, 'test_peaks.bed'), format='bed-summit')
