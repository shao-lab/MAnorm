Welcome to MAnorm
=================

.. image:: https://travis-ci.org/shao-lab/MAnorm.svg?branch=master
   :alt: Travis Build
   :target: https://travis-ci.org/shao-lab/MAnorm
.. image:: https://readthedocs.org/projects/manorm/badge/?version=latest
   :alt: Documentation Status
   :target: http://manorm.readthedocs.io/en/latest/?badge=latest
.. image:: https://img.shields.io/pypi/v/MAnorm.svg
   :alt: PyPI
   :target: https://pypi.python.org/pypi/MAnorm
.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
   :alt: Bioconda
   :target: http://bioconda.github.io/recipes/manorm/README.html
.. image:: https://img.shields.io/pypi/l/MAnorm.svg
   :alt: License
   :target: https://github.com/shao-lab/MAnorm/blob/master/LICENSE

What is MAnorm?
---------------

**MAnorm** is a robust model for quantitative comparison of ChIP-Seq data sets and you can use it for:

- Normalization of two ChIP-seq samples
- Quantitative comparison (differential analysis) of two ChIP-seq samples
- Evaluating the overlap enrichment of the protein binding sites(peaks)
- Elucidating underlying mechanisms of cell-type specific gene regulation

MAnorm is developed in Python with a command line interface, you can use it on Linux and macOS platforms.
Basically, it requires 4 input files (2 ChIP-seq peak files and 2 sequencing read files) corresponds to the
two samples under comparison. Starting from version 1.2.0, MAnorm supports multiple format for peaks/reads files.


Contents
--------

.. toctree::
   :maxdepth: 2

   intro
   tutorial
   changelog
   faq
   license
   contact

Citation
--------

If you use MAnorm or any derived code, please cite this paper in your publication:

`Shao Z, Zhang Y, Yuan GC, Orkin SH, Waxman DJ. (2012) MAnorm: a robust model for quantitative comparison of ChIP-Seq data sets. Genome Biol. Mar 16;13(3):R16.`__

.. __: http://genomebiology.com/2012/13/3/R16/abstract


---------------

The Python version of MAnorm is developed by ShaoLab_ at `CAS-MPG Partner Institute for Computational Biology, SIBS, CAS`_.

.. seealso::
   GitHub repository of MAnorm: https://github.com/shao-lab/MAnorm

.. _ShaoLab: http://bioinfo.sibs.ac.cn/shaolab/
.. _CAS-MPG Partner Institute for Computational Biology, SIBS, CAS: http://www.picb.ac.cn/picb/indexeng.jsp

