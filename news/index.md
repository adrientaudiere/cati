# Changelog

## cati (development version)

## cati 0.99.6

CRAN release: 2026-03-31

## cati 0.99.5

- [`MinMaxMST()`](https://adrientaudiere.github.io/cati/reference/MinMaxMST.md)
  and
  [`SumBL()`](https://adrientaudiere.github.io/cati/reference/SumBL.md)
  now use
  [`cluster::daisy()`](https://rdrr.io/pkg/cluster/man/daisy.html)
  instead of `FD::gowdis()` for Gower distance computation, replacing
  the FD dependency which is scheduled for CRAN archival.
- [`plot.ComIndexMulti()`](https://adrientaudiere.github.io/cati/reference/ComIndexMulti.md),
  [`plot.listofindex()`](https://adrientaudiere.github.io/cati/reference/plot.listofindex.md),
  [`Pval()`](https://adrientaudiere.github.io/cati/reference/Pval.md),
  [`print.ComIndex()`](https://adrientaudiere.github.io/cati/reference/ComIndex.md),
  [`print.ComIndexMulti()`](https://adrientaudiere.github.io/cati/reference/ComIndexMulti.md),
  [`print.Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.md),
  [`ses()`](https://adrientaudiere.github.io/cati/reference/ses.md), and
  [`ses.listofindex()`](https://adrientaudiere.github.io/cati/reference/ses.listofindex.md)
  no longer error on R \>= 4.0 where
  [`class()`](https://rdrr.io/r/base/class.html) on a matrix returns two
  values (`"matrix" "array"`), and on R \>= 4.5 where such comparisons
  in `if` conditions are treated as errors.
- [`ses()`](https://adrientaudiere.github.io/cati/reference/ses.md) and
  [`ses.Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.md)
  no longer produce non-conformable array errors when the number of
  communities equals the number of traits (square observed matrix), as
  the permutation dimension of the null model is now detected
  automatically.
- [`sum_Tstats()`](https://adrientaudiere.github.io/cati/reference/Tstats.md)
  with `type = "percent"` or `type = "all"` no longer errors when traits
  contain `NA` values.
- Add test suite with testthat
