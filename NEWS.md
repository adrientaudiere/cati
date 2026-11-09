# cati 0.99.7

This release fixes the bugs found in an external audit of cati 0.99.6 (GH-11 to GH-21). Many of them change numerical results: the observed and null values of `Tstats()`, `ComIndex()` and `ComIndexMulti()`, and the SES and p-values derived from them. Results computed with earlier versions should be recomputed.

## Breaking changes

* `Pval()`: `pval.inf` and `pval.sup` are swapped to match `sum_Tstats()`, so `pval.inf` is now small when the observed value is lower than the null values. Scripts that read `pval.inf` or `pval.sup` silently get the other tail; swap the two names to keep the previous meaning (GH-19).
* `Pval()` now returns one element per trait instead of only the first trait (GH-20).
* `Tstats()`: rows of `T_IP.IC`, `T_IC.IR` and `T_PC.PR` follow `levels(ind.plot)`, the order already used by `reg.pool` lists, instead of the alphabetical order. When the levels of `ind.plot` are not alphabetical, the rows of `T_PC.PR` move; index rows by name rather than by position (GH-17).

## Bug fixes

* `barplot.decompCTRE()` works with `resume = FALSE` (GH-21).
* `ComIndex()` and `ComIndexMulti()` keep species, `by.factor` and traits aligned when the rows are not sorted by plot (GH-15).
* `ComIndex()` and `ComIndexMulti()` draw each trait from its own values in the `regional.pop` null model (GH-14).
* `ComIndex()` and `ComIndexMulti()` use the regional pool of each community when `reg.pool` is a list, apply `SE` in the `local` null model, accept `SE.reg.pool` as a list, report the number of species in `$sites_richness`, sample correctly in plots with a single individual, and give different draws to plots of equal size when `independantTraits = FALSE` (GH-21).
* `ComIndex()` and `ComIndexMulti()` work with a single trait.
* `ComIndexMulti()` builds the null values of each trait from that trait only (GH-13).
* `ComIndexMulti()` with `nullmodels = "regional.pop.prab"` no longer errors, and groups population-level values with `by.factor`.
* `CVNND()` no longer returns `NA` when `na.rm = FALSE` (GH-21).
* `decompCTRE()` excludes individuals with `NA` traits from the weights (GH-21).
* `MinMaxMST()`, `MinNND()`, `MNND()` and `SDNND()` count zero distances between tied trait values (GH-21).
* `plotDistri()` draws populations of a single individual at their value and builds its legend with any number of levels (GH-21).
* `plotRandtest()` reads the null values of `ComIndex` objects and of any `listofindex` in the right layout, and its per-trait plot now tests the mean of the index over communities (GH-21).
* `plotSESvar()` works with one or two traits and no longer changes the global palette (GH-21).
* `plotSpPop()` and `plotSpVar()` keep species and traits aligned when the rows are not sorted by plot (GH-15), match populations without splitting names on `_`, and no longer fail on `color.lm` with the default `min.ind.signif` (GH-21).
* `Pval()`: `pval.inf` and `pval.sup` are swapped to match `sum_Tstats()`, so `pval.inf` is now small when the observed value is lower than the null values. Scripts that read `pval.inf` or `pval.sup` silently get the other tail; swap the two names to keep the previous meaning (GH-19).
* `Pval()` now returns one element per trait instead of only the first trait (GH-20).
* `Pval()` and `sum_Tstats(type = "p.value")` count ties with the null values against significance, drop `NA` null values from the denominator, and return `NA` when the observed value or all null values are `NA` (GH-19).
* `RandCom()` gives distinct names to species beyond the 26th (GH-21).
* `RaoRel()` returns `FD$Beta_add` and `FD$Beta_prop` with `Jost = TRUE` (GH-21).
* `samplingSubsetData()` samples correctly from units of one individual with `type = "count"`, keeps the largest units first with `type = "factorBySize"`, and respects `MinSample` with `type = "propBySize"` (GH-21).
* `ses()` gains a `null.dims` argument to give the layout of a bare null array. `ses()`, `ses.listofindex()` and `Pval()` use the layout of `Tstats` and `ComIndex` null values instead of guessing it from sizes, which gave wrong values when the numbers of plots, traits or permutations were equal (GH-20).
* `sum_Tstats()`: the `T_IC.IR.inf` row of `type = "site"` and the `T_PC.PR.sup` row of `type = "binary"` (and of `type = "all"`) report their own statistic (GH-16).
* `Tstats()` matches the populations of a plot by exact name, so plot `p1` no longer pools the populations of `p10` and plot names may contain regex characters (GH-11).
* `Tstats()` null models 1, 2 and 2sp keep each permutation and plot together, so null values no longer mix individuals of other plots or other permutations (GH-12, GH-18).
* `Tstats()`: rows of `T_IP.IC`, `T_IC.IR` and `T_PC.PR` follow `levels(ind.plot)`, the order already used by `reg.pool` lists, instead of the alphabetical order. When the levels of `ind.plot` are not alphabetical, the rows of `T_PC.PR` move; index rows by name rather than by position (GH-17).
* `Tstats()` accepts a character `ind.plot`, unused levels, and `@` in species or plot names, and gives different draws to plots of equal size when `independantTraits = FALSE` (GH-21).
* `Tstats()` computes the population means of null model 2sp with `na.rm = TRUE`, and stores `NA` instead of the string `"NA"` for plots with a single individual (GH-12, GH-19).

# cati 0.99.6

# cati 0.99.5

* `MinMaxMST()` and `SumBL()` now use `cluster::daisy()` instead of `FD::gowdis()` for Gower distance computation, replacing the FD dependency which is scheduled for CRAN archival.
* `plot.ComIndexMulti()`, `plot.listofindex()`, `Pval()`, `print.ComIndex()`, `print.ComIndexMulti()`, `print.Tstats()`, `ses()`, and `ses.listofindex()` no longer error on R >= 4.0 where `class()` on a matrix returns two values (`"matrix" "array"`), and on R >= 4.5 where such comparisons in `if` conditions are treated as errors.
* `ses()` and `ses.Tstats()` no longer produce non-conformable array errors when the number of communities equals the number of traits (square observed matrix), as the permutation dimension of the null model is now detected automatically.
* `sum_Tstats()` with `type = "percent"` or `type = "all"` no longer errors when traits contain `NA` values.
* Add test suite with testthat
