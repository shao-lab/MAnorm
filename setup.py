import os
import re
import sys

from setuptools import find_packages, setup

py_version = sys.version_info[:2]
if py_version < (3, 6):
    raise RuntimeError("MAnorm requires Python 3.6+ to install!")

description = "A robust model for quantitative comparison of ChIP-Seq data " \
              "sets"

pkg_dir = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(pkg_dir, 'manorm', '__init__.py')) as fin:
    version = re.search(r"__version__ = '(.*?)'", fin.read()).group(1)

with open(os.path.join(pkg_dir, 'README.rst')) as fin:
    long_description = fin.read()


install_requires = [
    "numpy",
    "pysam>=0.15.0",
    "matplotlib>=3.0.0",
    "scikit-learn>=0.21.0",
]

extras_require = {
    "test": ["pytest>=4.0.0",
             "pytest-cov>=2.8.0"],
    "docs": ["sphinx>=2.0.0",
             "sphinx_rtd_theme"]
}

classifiers = [
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Operating System :: Unix",
    "Operating System :: POSIX",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering :: Bio-Informatics"
]

setup(
    name="MAnorm",
    version=version,
    description=description,
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url="https://github.com/shao-lab/MAnorm",
    project_urls={
        "Bug Tracker": "https://github.com/shao-lab/MAnorm/issues",
        "Documentation": "https://manorm.readthedocs.io",
        "Source Code": "https://github.com/shao-lab/MAnorm",
    },
    author="Hayden Sun",
    author_email="sunhongduo@picb.ac.cn",
    license="BSD",
    packages=find_packages(),
    entry_points={"console_scripts": ["manorm=manorm.cli:main"]},
    python_requires=">=3.6",
    install_requires=install_requires,
    extras_require=extras_require,
    classifiers=classifiers,
    zip_safe=False,
)
