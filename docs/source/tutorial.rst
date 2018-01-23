.. _tutorial:

========
Tutorial
========

.. contents::
   :local:

Installation
============

Like many other Python packages and bioinformatics softwares, MAnorm can be obtained easily from PyPI_ or Bioconda_.
The command below shows how to install the latest release of MAnorm in a convenient way, but you can also install it
from source code alternatively.

Prerequisites
-------------

.. tip::
   MAnorm is implemented under **Python 2.7** and will support **Python 3.X** in the following updates.

* **Python 2.7**
* setuptools
* numpy
* matplotlib
* statsmodels
* scipy

Install with pip
----------------
The latest release of MAnorm is available at PyPI_, you can install via ``pip``::

    $ pip install manorm

.. _PyPI: https://pypi.python.org/pypi/MAnorm

Install with conda
----------------------

You can also install MAnorm with conda_ through Bioconda_ channel::

   $ conda install -c bioconda manorm

.. _conda: https://conda.io/docs/
.. _Bioconda: https://bioconda.github.io/

Install from source code
------------------------

It's highly recommended to install MAnorm with ``pip`` or ``conda``. If you prefer to install it from source code,
please read the following steps:

The source code of MAnorm is hosted on GitHub_, and setuptools_ is required for installation.

.. _setuptools: https://setuptools.readthedocs.io/en/latest/
.. _GitHub: https://github.com/shao-lab/MAnorm

First, clone the repository of MAnorm::

   $ git clone https://github.com/shao-lab/MAnorm.git

Then, install MAnorm in the source directory::

   $ cd MAnorm
   $ python setup.py install

.. note::
   * You may need to install all dependencies listed in ``requirements.txt``.
   * You may need to modify ``$PATH`` and ``$PYTHONPATH`` manually to make it work.

Galaxy Installation
-------------------
MAnorm is available on Galaxy_, you can incorporate MAnorm into your own Galaxy instance.

Please search and install MAnorm via the `Galaxy Tool Shed`_.

.. _Galaxy: https://galaxyproject.org
.. _`Galaxy Tool Shed`: https://toolshed.g2.bx.psu.edu/view/haydensun/manorm

Usage of MAnorm
===============

To check whether MAnorm is properly installed, you can inspect the version of MAnorm by ``-v/--version`` option::

  $ manorm -v
  $ manorm --version

Command-Line Usage
------------------

MAnorm provide a console script ``manorm`` for running the program, the basic usage should as follows:

  $ manorm --p1 peaks_file1.xls --p2 peaks_file2.xls --r1 reads_file1.bed --r2 reads_file2.bed -o output_name

.. tip::
    Please use ``-h/--help`` for the details of all options.

Options
-------

-h, --help     Show help message and exit.
-v, --version  Show version number and exit.
--p1           **[Required]** Peaks file of sample1.
--p2           **[Required]** Peaks file of sample2.
--r1           **[Required]** Reads file of sample1.
--r2           **[Required]** Reads file of sample2.
--s1           Reads shiftsize of sample1. Default: 100
--s2           Reads shiftsize of sample2. Default: 100
-w             Width of window to calculate read density. Default: 1000
-d             Summit-to-summit distance cutoff for common peaks. Default: ``-w``/2
-n             Number of simulations to test the enrichment of peaks overlap between two samples.
-m             *M-value* cutoff to distinguish biased (sample-specific) peaks from unbiased peaks.
-p             *P-value* cutoff to define biased peaks.
-s             Output additional files which contains the results of original peaks.
--name1        Name of sample1. (experiment condition, cell-type etc.)
--name2        Name of sample2.
-o             **[Required]** Output directory.

**Further explanation:**

  * ``--s1/--s2``:
    These values are used to shift reads towards 3' direction to determine the precise binding site.
    Set as half of the fragment length.

  * ``-w``:
    Half of the window size when counting reads of the peak regions. MAnorm uses windows with unified length of
    2 * ``-w`` centered at peak summits/midpoints to calculate the read density. This value should match the typical
    length of peaks, a value of 1000 is recommended for sharp histone marks like H3K4me3 and H3K9/27ac, and 500 for
    transcription factors or DNase-Seq.

  * ``-d``:
    Summit-to-summit distance cutoff for common peaks. Default= ``-w`` / 2. Only overlapped peaks with summit-to-summit
    distance less than than this value are considered as real common peaks of two samples when fitting M-A normalization
    model.

  * ``-m``:
    `M-value` (log2 fold change) cutoff to distinguish biased peaks from unbiased peaks. Peaks with M-value >= ``-m``
    and P-value <= ``-p`` are defined as sample1-biased(specific) peaks, while peaks with M-value <= -1 * ``-m``
    and P-value <= ``-p`` are defined as sample2-biased peaks.

  * ``-s``:
    By default, MAnorm will write the comparison results of unique and merged common peaks in a single output file.
    With this option on, MAnorm will output two extra files which contains the results of the original(unmerged) peaks.

  * ``--name1/--name2``:
    If specified, it will be used to replace the peaks/reads input file name as the sample name in output files.

  * ``-o``:
    Output directory. When ``--name1`` and ``--name2`` are not specified, MAnorm will use it as the prefix of comparison
    output file.

Input Format
============

Format of Peaks file
--------------------

Standard **BED** format and **MACS xls** format are supported, other supported format are listed below::

  * 3-columns tab split format

    # chr   start end
      chr1  2345  4345
      chr1  3456  5456
      chr2  6543  8543

  * 4-columns tab split format

    # chr   start end   summit
      chr1  2345  4345  254
      chr1  3456  5456  127
      chr2  6543  8543  302

.. note::
   The fourth column **summit** is the relative position to **start**.


Format of Reads file
--------------------

Only **BED** format are supported for now. More format will be embedded in the following updates.


MAnorm Output
=============

1. output_name_all_MAvalues.xls

This is the main output result of MAnorm which contains the M-A values and normalized read density of each peak,
common peaks from two samples are merged together.

 * chr: chromosome name
 * start: start position of the peak
 * end: end position of the peak
 * summit: summit position of the peak (relative to start)
 * m_value: M value (log2 Fold change) of normalized read densities under comparison
 * a_value: A value (average signal strength) of normalized read densities under comparison
 * p_value
 * peak_group: indicates where the peak  is come from
 * normalized_read_density_in _sample1
 * normalized_read_density_in_sample2

 .. note::
    Coordinates in .xls file is under **1-based** coordinate-system.

2. output_filters/

  * sample1_biased_peaks.bed
  * sample2_biased_peaks.bed
  * output_name_unbiased_peaks.bed

3. output_tracks/

  * output_name_M_values.wig
  * output_name_A_values.wig
  * output_name_P_values.wig

4. output_figures/

  * output_name_MA_plot_before_normalization.png
  * output_name_MA_plot_after_normalization.png
  * output_name_MA_plot_with_P-value.png
  * output_name_read_density_on_common_peaks.png
