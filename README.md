# MAnorm
[![travis-ci](https://travis-ci.org/shao-lab/MAnorm.svg?branch=master)](https://travis-ci.org/shao-lab/MAnorm)
[![pypi](https://img.shields.io/pypi/v/MAnorm.svg)](https://pypi.python.org/pypi/MAnorm)
[![Documentation Status](https://readthedocs.org/projects/manorm/badge/?version=latest)](http://manorm.readthedocs.io/en/latest/?badge=latest)
[![license](https://img.shields.io/pypi/l/MAnorm.svg)](https://github.com/shao-lab/MAnorm/blob/master/LICENSE)

## Introduction

ChIP-Seq is widely used to characterize genome-wide binding patterns of transcription factors and other chromatin-associated proteins. Although comparison of ChIP-Seq data sets is critical for understanding cell type-dependent and cell state-specific binding, and thus the study of cell-specific gene regulation, few quantitative approaches have been developed.

Here, we present a simple and effective method, MAnorm, for quantitative comparison of ChIP-Seq data sets describing transcription factor binding sites and epigenetic modifications. The quantitative binding differences inferred by MAnorm showed strong correlation with both the changes in expression of target genes and the binding of cell type-specific regulators.

## Citation

> [MAnorm: a robust model for quantitative comparison of ChIP-Seq data sets.](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-3-r16) Shao Z, Zhang Y, Yuan GC, Orkin SH, Waxman DJ.   *Genome Biol. Mar 16;13(3):R16.*

## Installation

MAnorm uses [setuptools](https://setuptools.readthedocs.io/en/latest/) for installation from source code. The source code of MAnorm is hosted on GitHub: https://github.com/shao-lab/MAnorm

You can clone the repo and execute the following command under source directory: 
```
python setup.py install
```

The latest version release of MAnorm is also available at [PyPI](https://pypi.python.org/pypi/MAnorm): 
```
pip install manorm
```
Or you can install MAnorm via conda:

In preparation.

## Usage
manorm [-h] [-v] --p1 peaks_file1 --p2 peaks_file2 --r1 reads_file1 --r2 reads_file2 -o output_name

Example: 
```
manorm --p1 sample1_peaks.bed --p2 sample2_peaks.bed --r1 sample1_reads.bed --r2 sample2_reads.bed -o sample1_vs_sample2
```

**Note:** Using -h/--help for the details of all arguments.

## Peak Format
Standard BED format and MACS xls format are supported, other supported format are listed below.

### 3-columns tab split format
Chromosome Start End
```
chr1  2345  4345
chr1  3456  5456
chr2  6543  8543 
```

### 4-columns tab split format
Chromosome Start End Summit

**Note:** The fourth column **Summit** is the relative position to **Start**.
 
```
chr1  2345  4345  254
chr1  3456  5456  127
chr2  6543  8543  302
```

## Reads Format
Only BED format are supported for now. More format will be embedded in the following updates. 

## Output of MAnorm
1. Output_name_all_MAvalues.xls

This is the main output result of MAnorm which contains the M-A values and normalized read density of each peak.
The common peaks from two samples are merged together.
- chr: chromosome name
- start: start position of the peak
- end: end position of the peak
- summit: summit position of the peak (relative to start)
- m_value: M value (log2 Fold change) of normalized read densities under comparison  
- a_value: A value (average signal strength) of normalized read densities under comparison
- p_value 
- peak_group: indicates where the peak  is come from
- normalized_read_density_in _sample1
- normalized_read_density_in_sample2    

**Note**: Coordinates in .xls file is under **1-based** coordinate-system.

2. output_filters

- sample1_biased_peaks.bed
- sample2_biased_peaks.bed
- output_name_unbiased_peaks.bed

3. output_tracks
- output_name_M_values.wig
- output_name_A_values.wig
- output_name_P_values.wig

4. output_figures
- output_name_MA_plot_before_normalization.png
- output_name_MA_plot_after_normalization.png
- output_name_MA_plot_with_P-value.png
- output_name_read_density_on_common_peaks.png
  
## License

[BSD 3-Clause License](https://github.com/shao-lab/MAnorm/blob/master/LICENSE)
