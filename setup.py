#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import re
import sys

from setuptools import find_packages, setup

py_version = sys.version_info[:2]
if not (py_version == (2, 7) or py_version >= (3, 4)):
    raise RuntimeError("MAnorm requires Python 2.7 or 3.4+ to install!")

description = "A robust model for quantitative comparison of ChIP-Seq data sets."

pkg_dir = os.path.abspath(os.path.dirname(__file__))

with io.open(os.path.join(pkg_dir, 'manorm', '__init__.py'), encoding='utf-8') as fin:
    version = re.search(r'__version__ = \'(.*?)\'', fin.read()).group(1)

with io.open(os.path.join(pkg_dir, 'README.rst'), encoding='utf-8') as fin:
    long_description = fin.read()


install_requires = [
    "setuptools",
    "numpy",
    "pysam>=0.15.0",
    "matplotlib>=2.2.3",
    "scikit-learn>=0.19.2",
]

extras_require = {
    "test": ["pytest >= 3.0.0",
             "pytest-cov >= 2.6.0"],
    "docs": ["sphinx >= 1.8.0",
             "sphinx_rtd_theme"]
}

classifiers = [
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Operating System :: POSIX",
    "Operating System :: MacOS :: MacOS X"
    "Programming Language :: Python",
    "Programming Language :: Python :: 2",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Topic :: Scientific/Engineering :: Bio-Informatics"
]

setup(
    name="MAnorm",
    version=version,
    description=description,
    long_description=long_description,
    author="Hayden Sun",
    author_email="sunhongduo@picb.ac.cn",
    url="https://github.com/shao-lab/MAnorm",
    license="BSD",
    packages=find_packages(),
    entry_points={"console_scripts": ["manorm=manorm.cli:main"]},
    python_requires=">=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*",
    install_requires=install_requires,
    extras_require=extras_require,
    classifiers=classifiers,
    zip_safe=False,
)
