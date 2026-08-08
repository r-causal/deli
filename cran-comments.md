## Resubmission

This is a resubmission. The previous submission was returned with the comment
that the package has examples for unexported functions. In this version:

* The methods documented on the estimator method topic pages (`coef()`,
  `vcov()`, `summary()`, `predict()`, `tidy()`, `augment()`, and the other
  standard accessors) are methods for an S7 class, and they were registered
  only at package load time by `S7::methods_register()`. Nothing about them
  appeared in NAMESPACE, so although the documentation and its examples were
  accurate, the methods read as unexported.
* Each of those methods now has its own `S3method()` directive in NAMESPACE,
  one per method. Every function documented with examples is therefore either
  exported or formally registered as a method.
* The examples themselves are unchanged.

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
