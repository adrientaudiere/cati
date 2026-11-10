# cati 0.99.7

This release fixes the bugs found in an external audit of cati 0.99.6 (GH-11 to GH-21). Many of them change numerical results. The section below says which analyses are affected and by how much; the figures come from `finch.ind` (6 plots, 4 traits, 199 permutations), comparing 0.99.6 and 0.99.7. Between two seeds of 0.99.7, the SES of `Tstats()` vary by 0.5 (median), so smaller differences are noise.

## Which results change

Results that do not change:

* Observed values of `Tstats()`, when no plot name is contained in another one and the levels of `ind.plot` are in alphabetical order.
* `ComIndex()` with the `local` and `regional.ind` null models, when the rows are sorted by plot, `reg.pool` is not a list, `SE = 0`, `independantTraits = TRUE` and the index does not compute nearest-neighbour distances on tied values. On `finch.ind`, observed values, null values and SES are identical.
* Null values of `T_IC.IR` and `T_PC.PR` in `Tstats()` when `reg.pool` is not a list, the rows are sorted by plot and the traits have no `NA`: the differences are within noise.
* P-values of `sum_Tstats()`, except when the observed value equals null values or is `NA`, or when some null values are `NA`.

Results that change:

* `Tstats()`, null values of `T_IP.IC`, in every analysis (GH-12). On `finch.ind`, the SES change by 15% in median (up to 77%) and keep their sign, and all 23 significant cells stay significant. The change is larger when the rows are not sorted by plot.
* `Tstats()`, null values of `T_PC.PR` when the traits contain `NA` (GH-19): populations with at least one `NA` value were left out of the regional pool. On `finch.ind` (1137 `NA` values), 7 of the 24 cells become significant (SES 2 to 3.6), against none before.
* `Tstats()`, null values of `T_IC.IR` when `reg.pool` is a list (GH-18): each community was compared with a mix of all the regional pools.
* `Tstats()`, observed `T_IP.IC` when a plot name is contained in another one (`1` and `10`, `p1` and `p10`) or contains regex characters (GH-11): the populations of several plots were pooled (0.44 instead of 0.006 in the GH-11 example).
* `Tstats()`, rows of `T_IP.IC` and `T_IC.IR` when the levels of `ind.plot` are not in alphabetical order, for example a numeric plot code (GH-17): values were labelled with the name of another plot, so any comparison with plot covariates used the wrong plots, and `sum_Tstats()` and `barplot.Tstats()` combined statistics of different plots.
* `ComIndex()` and `ComIndexMulti()` with the `regional.pop` null model and more than one trait (GH-14): the null values were drawn from all traits pooled. On `finch.ind`, the significance changes in 23 of 24 cells for the mean and 19 of 24 for the variance.
* `ComIndexMulti()` with more than one trait, with any null model (GH-13): the null values mixed the traits. On `finch.ind`, the `local` null values of the mean of the first trait ranged from 9.6 to 70.7 instead of 67.0 to 68.3.
* `ComIndex()` and `ComIndexMulti()` when the rows are not sorted by plot (GH-15): indices using `sp` or `by.factor`, and the observed values of `regional.pop`, were computed on misaligned species.
* `MNND()`, `MinNND()`, `SDNND()` and `MinMaxMST()` when trait values are tied (GH-21): zero distances were dropped. Tied values are the rule in `finch.ind` (92% to 98% of the values), and the significance of `MNND` under `regional.ind` changes in 13 of 24 cells.
* `ComIndex()` and `ComIndexMulti()` null values when `reg.pool` is a list, `SE` is not 0 with the `local` null model, or `independantTraits = FALSE` (GH-21).
* `ses()`, `ses.Tstats()`, `plot.Tstats()` and `Pval()` when the number of plots equals the number of traits, or, for `ComIndex()`, when the number of permutations equals one of them (GH-20): rows and columns were swapped. For a `ComIndex()` index that returns one value per trait, such as `"mean(x)"`, the layout is still guessed from sizes.
* `Pval()`: `pval.inf` and `pval.sup` are swapped, and the traits after the first one were missing (GH-19, GH-20).
* `sum_Tstats()`: the site list of `T_IC.IR.inf` and the binary row of `T_PC.PR.sup` showed another statistic (GH-16).

## Breaking changes

* `Pval()`: `pval.inf` and `pval.sup` are swapped to match `sum_Tstats()`, so `pval.inf` is now small when the observed value is lower than the null values. Scripts that read `pval.inf` or `pval.sup` silently get the other tail; swap the two names to keep the previous meaning (GH-19).
* `Pval()` now returns one element per trait instead of only the first trait (GH-20).
* `Tstats()`: rows of `T_IP.IC`, `T_IC.IR` and `T_PC.PR` follow `levels(ind.plot)`, the order already used by `reg.pool` lists, instead of the alphabetical order. When the levels of `ind.plot` are not alphabetical, the rows of `T_PC.PR` move; index rows by name rather than by position (GH-17).
* `MNND()`, `MinNND()`, `SDNND()` and `MinMaxMST()` keep the zero distances between individuals with the same trait values, and only exclude the distance of an individual to itself, as `CVNND()` already did. They used to drop every zero distance, which ignored tied values. With traits measured at a coarse resolution, `MinNND()` and `MinMaxMST()` are now often 0 (GH-21).

## Bug fixes

* `barplot.decompCTRE()` works with `resume = FALSE` (GH-21).
* `ComIndex()` and `ComIndexMulti()` keep species, `by.factor` and traits aligned when the rows are not sorted by plot (GH-15).
* `ComIndex()` and `ComIndexMulti()` draw each trait from its own values in the `regional.pop` null model (GH-14).
* `ComIndex()` and `ComIndexMulti()` use the regional pool of each community when `reg.pool` is a list, apply `SE` in the `local` null model, accept `SE.reg.pool` as a list, report the number of species in `$sites_richness`, sample correctly in plots with a single individual, and give different draws to plots of equal size when `independantTraits = FALSE` (GH-21).
* `ComIndex()` and `ComIndexMulti()` work with a single trait.
* `ComIndex()` with `nullmodels = "regional.pop"` or `"regional.pop.prab"` evaluates `index` on one value per population with `ind.plot` and `sp` giving the plot and species of each population, so indices such as `"tapply(x, ind.plot, mean)"` work with every null model. The `ComIndex()` example no longer needs a workspace variable.
* `ComIndexMulti()` builds the null values of each trait from that trait only (GH-13).
* `ComIndexMulti()` with `nullmodels = "regional.pop.prab"` no longer errors, and groups population-level values with `by.factor`.
* `CVNND()` no longer returns `NA` when `na.rm = FALSE` (GH-21).
* `decompCTRE()` excludes individuals with `NA` traits from the weights (GH-21).
* `plotDistri()` draws populations of a single individual at their value and builds its legend with any number of levels (GH-21).
* `plotRandtest()` reads the null values of `Tstats` objects and of `ComIndex` indices that return one value per community in the right layout, and its per-trait plot now tests the mean of the index over communities (GH-21). It still fails on `ComIndexMulti` objects and on `ComIndex` indices that return one value per trait.
* `plotSESvar()` works with one or two traits and no longer changes the global palette (GH-21).
* `plotSpPop()` and `plotSpVar()` keep species and traits aligned when the rows are not sorted by plot (GH-15), match populations without splitting names on `_`, and no longer fail on `color.lm` with the default `min.ind.signif` (GH-21).
* `Pval()` and `sum_Tstats(type = "p.value")` count ties with the null values against significance, drop `NA` null values from the denominator, and return `NA` when the observed value or all null values are `NA` (GH-19).
* `RandCom()` gives distinct names to species beyond the 26th (GH-21).
* `RaoRel()` returns `FD$Beta_add` and `FD$Beta_prop` with `Jost = TRUE` (GH-21).
* `samplingSubsetData()` samples correctly from units of one individual with `type = "count"`, keeps the largest units first with `type = "factorBySize"`, and respects `MinSample` with `type = "propBySize"` (GH-21).
* `ses()` gains a `null.dims` argument to give the layout of a bare null array. `ses()`, `ses.listofindex()` and `Pval()` use the layout of `Tstats` and `ComIndex` null values instead of guessing it from sizes, which gave wrong values when the numbers of plots, traits or permutations were equal (GH-20). The layout is still guessed for `ComIndex` indices that return one value per trait.
* `sum_Tstats()`: the `T_IC.IR.inf` row of `type = "site"` and the `T_PC.PR.sup` row of `type = "binary"` (and of `type = "all"`) report their own statistic (GH-16).
* `Tstats()` matches the populations of a plot by exact name, so plot `p1` no longer pools the populations of `p10` and plot names may contain regex characters (GH-11).
* `Tstats()` null models 1, 2 and 2sp keep each permutation and plot together, so null values no longer mix individuals of other plots or other permutations (GH-12, GH-18).
* `Tstats()` accepts a character `ind.plot`, unused levels, and `@` in species or plot names, and gives different draws to plots of equal size when `independantTraits = FALSE` (GH-21).
* `Tstats()` computes the population means of null model 2sp with `na.rm = TRUE`, and stores `NA` instead of the string `"NA"` for plots with a single individual (GH-12, GH-19).

# cati 0.99.6

# cati 0.99.5

* `MinMaxMST()` and `SumBL()` now use `cluster::daisy()` instead of `FD::gowdis()` for Gower distance computation, replacing the FD dependency which is scheduled for CRAN archival.
* `plot.ComIndexMulti()`, `plot.listofindex()`, `Pval()`, `print.ComIndex()`, `print.ComIndexMulti()`, `print.Tstats()`, `ses()`, and `ses.listofindex()` no longer error on R >= 4.0 where `class()` on a matrix returns two values (`"matrix" "array"`), and on R >= 4.5 where such comparisons in `if` conditions are treated as errors.
* `ses()` and `ses.Tstats()` no longer produce non-conformable array errors when the number of communities equals the number of traits (square observed matrix), as the permutation dimension of the null model is now detected automatically.
* `sum_Tstats()` with `type = "percent"` or `type = "all"` no longer errors when traits contain `NA` values.
* Add test suite with testthat
