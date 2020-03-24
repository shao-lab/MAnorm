import os

import pytest

from manorm.exceptions import FileFormatError
from manorm.region import load_manorm_peaks


def test_unsupported_format(data_dir):
    with pytest.raises(ValueError):
        load_manorm_peaks(os.path.join(data_dir, 'test_peaks.bed'),
                          format='unknown_format')


def test_invalid_format(data_dir):
    with pytest.raises(FileFormatError):
        load_manorm_peaks(os.path.join(data_dir, 'test_peaks.bed'),
                          format='bed3-summit')
