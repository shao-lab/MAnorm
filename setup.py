#!/usr/bin/env python

from setuptools import setup, find_packages
from manorm import __version__

DESCRIPTION = "A robust model for quantitative comparison of ChIP-Seq data sets."
with open('README.rst') as fin:
    LONG_DESCRIPTION = fin.read()

INSTALL_REQUIRES = ["numpy",
                    "matplotlib",
                    "statsmodels"]

CLASSIFIERS = ["Development Status :: 5 - Production/Stable",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: BSD License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Programming Language :: Python :: 2",
               "Programming Language :: Python :: 2.7",
               "Topic :: Scientific/Engineering :: Bio-Informatics"]

setup(
    name="MAnorm",
    version=__version__,
    packages=find_packages(exclude=["test"]),
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    url="https://github.com/shao-lab/MAnorm",
    author="Semal",
    author_email="",
    maintainer="Hayden Sun",
    maintainer_email="sunhongduo@picb.ac.cn",
    license="BSD",
    entry_points={"console_scripts": ["manorm=manorm.cmdline:main"]},
    install_requires=INSTALL_REQUIRES,
    classifiers=CLASSIFIERS,
    zip_safe=False,
)
