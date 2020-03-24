import os

import pytest

from manorm.exceptions import FileFormatError, FormatModeConflictError
from manorm.read import load_reads


def test_matched_bed(data_dir):
    with pytest.raises(FormatModeConflictError):
        load_reads(os.path.join(data_dir, 'test_reads.bed'), format='bed',
                   paired=True)


def test_matched_bedpe(data_dir):
    with pytest.raises(FormatModeConflictError):
        load_reads(os.path.join(data_dir, 'test_reads.bedpe'), format='bedpe',
                   paired=False)


def test_unsupported_format(data_dir):
    with pytest.raises(ValueError):
        load_reads(os.path.join(data_dir, 'test_reads.bed'),
                   format='unknown_format')


def test_invalid_format(data_dir):
    with pytest.raises(FileFormatError):
        load_reads(os.path.join(data_dir, 'test_reads.bed'), format='bedpe',
                   paired=True)
