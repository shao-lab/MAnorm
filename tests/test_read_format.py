import os

import pytest

from manorm.exceptions import InvalidFormatError, UnmatchedBedFormatError, UnsupportedFormatError
from manorm.read import load_reads


def test_matched_bed(data_dir):
    with pytest.raises(UnmatchedBedFormatError):
        reads = load_reads(os.path.join(data_dir, 'test_reads.bed'), format='bed', paired=True)


def test_matched_bedpe(data_dir):
    with pytest.raises(UnmatchedBedFormatError):
        reads = load_reads(os.path.join(data_dir, 'test_reads.bedpe'), format='bedpe', paired=False)


def test_unsupported_format(data_dir):
    with pytest.raises(UnsupportedFormatError):
        reads = load_reads(os.path.join(data_dir, 'test_reads.bed'), format='unknown_format')


def test_invalid_format(data_dir):
    with pytest.raises(InvalidFormatError):
        reads = load_reads(os.path.join(data_dir, 'test_reads.bed'), format='bedpe', paired=True)