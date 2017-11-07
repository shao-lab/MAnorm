#!/usr/bin/env python

from setuptools import setup, find_packages
from manorm import __version__

DESCRIPTION = "A robust model for quantitative comparison of ChIP-Seq data sets."
LONG_DESCRIPTION = """ChIP-Seq is widely used to characterize genome-wide binding patterns of transcription factors 
and other chromatin-associated proteins. Although comparison of ChIP-Seq data sets is critical for understanding the 
role of their cell type-specific binding on modulating gene regulation programs, few quantitative approaches have 
been developed. Here, we present a simple and effective method, MAnorm, for quantitative comparison of ChIP-Seq data 
sets describing transcription factor binding sites and epigenetic modifications. The quantitative binding differences 
inferred by MAnorm showed strong correlation with both the changes in expression of target genes and the binding of 
cell type-specific regulators. 

MAnorm uses common peaks of two samples as a reference to build the rescaling model for normalization, which is based 
on the empirical assumption that if a chromatin-associated protein has a substantial number of peaks shared in two 
conditions, the binding at these common regions will tend to be determined by similar mechanisms, and thus should  
exhibit similar global binding intensities across samples. 

The normalized M value given by MAnorm was used as a quantitative measure of differential binding in each peak region 
between two samples, with peak regions associated with larger absolute M values exhibiting greater binding differences 
between two samples.  

MAnorm exhibited excellent performance in quantitative comparison of ChIP-Seq data sets for both epigenetic 
modifications and transcription factors. The quantitative binding differences inferred by MAnorm were highly 
correlated with both the changes in expression of target genes and also the binding of cell type-specific regulators. 
With the accumulation of ChIP-seq data sets, MAnorm should serve as a powerful tool for obtaining a more 
comprehensive understanding of cell type-specific and cell state-specific regulation during organism development and 
disease onset. 

**Citation** MAnorm: a robust model for quantitative comparison of ChIP-Seq data sets. Shao Z, Zhang Y, Yuan GC, Orkin SH, Waxman DJ. Genome Biol. Mar 16;13(3):R16. 


"""

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
