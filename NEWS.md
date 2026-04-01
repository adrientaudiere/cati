# cati (development version)

# cati 0.99.6

# cati 0.99.5

* `MinMaxMST()` and `SumBL()` now use `cluster::daisy()` instead of `FD::gowdis()` for Gower distance computation, replacing the FD dependency which is scheduled for CRAN archival.
* `plot.ComIndexMulti()`, `plot.listofindex()`, `Pval()`, `print.ComIndex()`, `print.ComIndexMulti()`, `print.Tstats()`, `ses()`, and `ses.listofindex()` no longer error on R >= 4.0 where `class()` on a matrix returns two values (`"matrix" "array"`), and on R >= 4.5 where such comparisons in `if` conditions are treated as errors.
* `ses()` and `ses.Tstats()` no longer produce non-conformable array errors when the number of communities equals the number of traits (square observed matrix), as the permutation dimension of the null model is now detected automatically.
* `sum_Tstats()` with `type = "percent"` or `type = "all"` no longer errors when traits contain `NA` values.
* Add test suite with testthat
