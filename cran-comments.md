## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Data and test file provenance

* The datasets in `data/` are re-transcribed or imported from the MIT-licensed
  Python 'delicatessen' library and from published sources. Each one is
  documented with references in its .Rd file.
* `tests/testthat/fixtures/` holds about 2.7 MB of JSON reference output
  generated with Python 'delicatessen', along with the Python scripts used to
  generate most of it. The fixtures cross-validate this package's numerical
  results against the original library.
* Paul Zivich, the author of 'delicatessen', is credited in Authors@R as a
  contributor, and the library is cited in the Description field.

## References

Method references are given in the Description field with DOIs.
