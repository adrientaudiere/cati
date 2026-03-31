## Resubmission

This is a resubmission. Changes in this version:

* Replaced the `FD` dependency (scheduled for CRAN archival on 2026-04-15)
  with `cluster::daisy()` in `SumBL()` and `MinMaxMST()`.
* Fixed R >= 4.5 compatibility: `class(matrix) == "string"` comparisons in
  `if` conditions now cause errors in R 4.5; replaced with `is.list()`,
  `inherits()`, and `identical(class(x), "list")` throughout.
* Fixed `ses()` / `ses.Tstats()`: non-conformable array error when the number
  of communities equals the number of traits.
* Fixed `sum_Tstats()` type "percent"/"all": NA subscript error with missing
  trait values.

The maintainer email has changed to adrien.taudiere@zaclys.net.

## R CMD check results

0 errors | 0 warnings | 2 notes

* NOTE: New maintainer email address (intentional change).
* NOTE: "unable to verify current time" — network issue in local check
  environment; not reproducible on CRAN infrastructure.

## Downstream dependencies

There are no downstream dependencies for this package.
