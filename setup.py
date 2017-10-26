from setuptools import setup, find_packages
from manorm import __version__

INSTALL_REQUIRES = ["numpy",
                    "scipy",
                    "matplotlib",
                    "statsmodels"
                    ]

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Education",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
               "Operating System :: POSIX :: Linux",
               "Operating System :: Unix",
               "Programming Language :: Python",
               "Programming Language :: Python :: 2",
               "Programming Language :: Python :: 2.7",
               "Topic :: Scientific/Engineering :: Bio-Informatics",
               ]

with open('README.md', 'r') as fin:
    long_description = fin.read()

setup(
    name="MAnorm",
    version=__version__,
    packages=find_packages(),
    scripts={'bin/manorm'},
    description="MAnorm -- A robust model for quantitative comparison of ChIP-Seq data sets, developed by Shao lab.",
    long_description=long_description,
    url="http://bioinfo.sibs.ac.cn/shaolab/MAnorm/MAnorm.htm",
    author="Semal",
    author_email="gzhsss2@gmail.com",
    maintainer="Hayden Sun",
    maintainer_email="sunhongduo@picb.ac.cn",
    license="GPLv3",
    install_requires=INSTALL_REQUIRES,
    classifiers=CLASSIFIERS,
    zip_safe=False,
    )
