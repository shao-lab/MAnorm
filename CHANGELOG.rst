ChangeLog
=========

v1.3.0 (2020-xx-xx)
-------------------

* Drop support for Python versions under 3.6
* Support broadPeak format
* Refactor for better compatbility with MotifScan and MAmotif

v1.2.0 (2019-01-03)
-------------------

* Support python3.4+
* Support MACS2-xls/narrowPeak format
* Support BEDPE/SAM/BAM format and paired-end mode for sequencing reads
* Drop support for Windows platform
* Replace statsmodels with scikit-learn for robust linear regression
* Fix a bug in summit calculation
* Various performance improvements


v1.1.4 (2018-08-17)
-------------------

* Fix an issue in setting matplotlib backend


v1.1.3 (2018-01-19)
-------------------

* Fix a bug in the file name of filtered biased peaks

* Fix a typo


v1.1.2 (2018-01-18)
-------------------

* Keep five digits for floats in the output files

* Fix a typo


v1.1.1 (2018-01-11)
-------------------

* Add test module

Bugs fixed:

* Rename the file name of filtered biased peaks


v1.1 (2017-11-07)
-----------------

Improvements:

* Refactor the package for better performance and compatibility

Bugs fixed:

* Fix the coordinates of peaks to be consistent with the corresponding coordinate system
* Fix the approximate equation in p-value calculation
* Fix the summit calculation of merged common peaks
