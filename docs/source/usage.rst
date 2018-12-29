.. _usage:

=====
Usage
=====

.. contents::
   :local:

Usage of MAnorm
===============

Command Line Interface
----------------------

First, to check whether MAnorm is properly installed, you can inspect the version of MAnorm by ``-v/--version`` option:

.. code-block:: shell

  $ manorm --version

Basic Usage
-----------

MAnorm provides a console script ``manorm`` for running the program, the basic usage is as follows:

.. code-block:: shell

  $ manorm --p1 peaks_file1.bed --p2 peaks_file2.bed --r1 reads_file1.bed --r2 reads_file2.bed
  --n1 name1 --n2 name2 --dir output_dir

.. tip::
    Please use ``-h/--help`` for the details of all options.

Options
-------

-h, --help           Show help message and exit.
-v, --version        Show version number and exit.
--p1, --peak1        **[Required]** Peak file of sample1.
--p2, --peak2        **[Required]** Peak file of sample2.
--pf, --peak-format  Format of peak files. Default: bed
--r1, --read1        **[Required]** Read file of sample1.
--r2, --read2        **[Required]** Read file of sample2.
--rf, --read-format  Format of read files. Default: bed
--n1, --name1        Name of sample 1.
--n2, --name2        Name of sample 2.
--s1, --shiftsize1   Reads shiftsize of sample1. Default: 100
--s2, --shiftsize2   Reads shiftsize of sample2. Default: 100
--pe, --paired-end   Paired-end mode.
-w, --window-size    Window size to count reads and calculate read density. Default: 2000
--summit-dis         Summit-to-summit distance  cutoff for common peaks. Default: ``-w``/4
--n-random           Number of simulations to test the enrichment of peaks overlap between two samples.
-m, --m-cutoff       Absolute *M* value (*log*:sub:`2`-ratio) cutoff to define biased (differential binding) peaks.
-p, --p-cutoff       *P* value cutoff to define biased peaks.
--oa, --output-all   Output additional files which contains the results of original (unmerged) peaks.
--dir                **[Required]** Output directory.
--verbose            Enable verbose log messages.

**Detailed explanation on options:**

  * ``--pf/--peak-format``:

    The format of peak files. For more details, see :ref:`peak file formats`.

  * ``--rf/--read-format``:

    The format of read files. For more details, see :ref:`read file formats`.

  * ``--n1/--name1`` and ``--n2/--name2``:

    These two options specify the sample names which are used in all output files.
    If not specified, the name of the peak file will be used as the sample name.

  * ``--s1/--shiftsize1`` and ``--s2/--shiftsize2``:

    These values are used to shift **single-end** reads towards 3' direction and the 5' end of each
    shifted read is used to represent the genomic locus of underlying DNA fragment. Set to half
    of DNA fragment size of the ChIP-seq library. These options are disabled in paired-end mode.

  * ``--pe/--paired-end``:

    Paired-end mode. The middle point of each read pair is used to represent the genomic locus of
    underlying DNA fragment. ``--s1`` and ``--s2`` are ignored with this option on.

  * ``-w/--window-size``:

    Window size to count reads and calculate read densities. 2000 is recommended for sharp histone
    marks like H3K4me3 and H3K27ac, and 1000 for TFs or DNase-seq. Default: 2000

  * ``--summit-dis``:

    Overlapping common peaks with summit-to-summit distance beyond this are excluded in model fitting.
    This option is used to exclude common peaks that only overlap on the edge of each other.
    Default: ``-w/--window-size``/4

  * ``--oa/--output-all``:

    By default, MAnorm only write the comparison results of unique and merged common peaks in a single
    output file. With this option on, MAnorm will write two extra files which contains the results of
    the original(unmerged) peaks.


Input File Format
=================

.. _`peak file formats`:

Format of peak files
--------------------

BED format
^^^^^^^^^^

Standard `BED`_ format is supported, the first 3 columns (``chrom``, ``start``, ``end``) of the bed file are used.


MACS
^^^^

`MACS`_ xls format is supported, MAnorm uses the ``chrom``, ``start``, ``end`` and ``summit`` information.

MACS2
^^^^^

`MACS2`_ xls format is supported, MAnorm uses the ``chrom``, ``start``, ``end`` and ``summit`` information.

narrowPeak
^^^^^^^^^^

ENCODE `narrowPeak`_ format is supported, the first 3 columns (``chrom``, ``start``, ``end``) are used and if
the 10th column is available, MAnorm uses it as the ``summit`` coordinate.

BED-summit
^^^^^^^^^^

A customized BED format named as ``BED-summit`` is also supported, the first 3 columns is same as ``BED`` format,
but the 4th columns should be the ``summit`` position (**relative** position to ``start``)

.. _`read file formats`:

Format of read files
--------------------

.. note:: MAnorm does not excluded any duplicated reads, and you may need use other tools to remove
          duplicates in advance to if you want.

BED format
^^^^^^^^^^

.. note:: BED format can only be used in **single-end** mode.

Standard `BED`_ format is supported.

BEDPE format
^^^^^^^^^^^^

.. note:: BEDPE format can only be used in **paired-end** mode.

`BEDPE`_ format which is defined by `bedtools`_ is also supported. Paired reads with
both ends mapped to a same chromosome are counted.

SAM/BAM format
^^^^^^^^^^^^^^

Standard `SAM`_ format and its binary form `BAM`_ format are supported.

When in paired-end mode, only proper paired mapped reads with both ends mapped to
the same chromosome are counted.

MAnorm Output Files
===================

1. *_all_MAvalues.xls

This is the main output result of MAnorm which contains the M-A values and normalized
read density of each peak, common peaks from two samples are merged together.

 - chr: chromosome name
 - start: start position of the peak
 - end: end position of the peak
 - summit: summit position of the peak (relative to start)
 - m_value: *M* value (*log*:sub:`2` fold change) of normalized read densities under comparison
 - a_value: *A* value (average signal strength) of normalized read densities under comparison
 - p_value
 - peak_group: indicates where the peak is come from and whether it is a common peak
 - normalized_read_density_in_sample1
 - normalized_read_density_in_sample2

 .. note::
    Coordinates in .xls file is under **1-based** coordinate-system.

2. output_filters/

This folder contains the filtered biased/unbiased peaks in BED format.

  - \*_M_above_*_biased_peaks.bed
  - \*_M_below_*_biased_peaks.bed
  - \*_unbiased_peaks.bed

3. output_tracks/

These files are genome track files of M values, A values and P values in ``wig`` format,
you can upload these files to a genome browser to visualize them.

  - \*_M_values.wig
  - \*_A_values.wig
  - \*_P_values.wig

4. output_figures/

This folder contains M-A plots before/after normalization and a scatter plot which shows the
scaling relationship between two samples.

  - \*_MA_plot_before_normalization.png
  - \*_MA_plot_after_normalization.png
  - \*_MA_plot_with_P_value.png
  - \*_read_density_on_common_peaks.png

.. _BED: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
.. _MACS: http://liulab.dfci.harvard.edu/MACS/README.html
.. _MACS2: https://github.com/taoliu/MACS
.. _narrowPeak: https://genome.ucsc.edu/FAQ/FAQformat.html#format12
.. _BEDPE: https://bedtools.readthedocs.io/en/latest/content/general-usage.html
.. _SAM: https://samtools.github.io/hts-specs/SAMv1.pdf
.. _BAM: https://samtools.github.io/hts-specs/SAMv1.pdf
.. _bedtools: https://bedtools.readthedocs.io/en/latest/index.html
