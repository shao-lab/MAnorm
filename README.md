# MAnorm
## Introduction

ChIP-Seq is widely used to characterize genome-wide binding patterns of transcription factors and other chromatin-associated proteins. Although comparison of ChIP-Seq data sets is critical for understanding cell type-dependent and cell state-specific binding, and thus the study of cell-specific gene regulation, few quantitative approaches have been developed.

Here, we present a simple and effective method, MAnorm, for quantitative comparison of ChIP-Seq data sets describing transcription factor binding sites and epigenetic modifications. The quantitative binding differences inferred by MAnorm showed strong correlation with both the changes in expression of target genes and the binding of cell type-specific regulators.

## Citation

> [MAnorm: a robust model for quantitative comparison of ChIP-Seq data sets.](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2012-13-3-r16) Shao Z, Zhang Y, Yuan GC, Orkin SH, Waxman DJ.   *Genome Biol. Mar 16;13(3):R16.*

## Installation
### pip
### conda
### setup.py

## Command
manorm --p1 peak_file1 --p2 peak_file2 --r1 reads_file1 --r2 reads_file2 -o output_name

**Note:** Using -h/--help for details.

## Peak Format
Standard BED format and MACS xls format are supported. 
Other supported format are listed below.

### 3-columns tab split format
Chromosome Start End
```
chr1  2345  4345
chr1  3456  5456
chr2  6543  8543 
```

### 4-columns tab split format
Chromosome Start End Summit

**Note:** The fourth column 'summit' is the relative position to 'start'.
 
```
chr1  2345  4345  254
chr1  3456  5456  127
chr2  6543  8543  302
```

## Reads Format
Only BED format are supported for now. More format will be embedded in the following updates. 

## License

## Contact
