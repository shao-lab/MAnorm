import os
import shutil
import pytest


@pytest.fixture(scope='session')
def data_dir():
    """Return the directory of data files."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')


@pytest.fixture(scope='function')
def tmp_dir(data_dir):
    """Generate a temporary output directory to write output files."""
    tmp_dir = os.path.join(data_dir, 'manorm_tmp_output')
    yield tmp_dir
    shutil.rmtree(tmp_dir)
