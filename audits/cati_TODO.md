# Audit Report: cati v0.99.5

Date: 2026-03-30

## Summary

- Total exported functions: ~40 (grouped across source file sections)
- Tested: 80 test cases
- Skipped: 0
- OK: 73
- Warnings: 0
- Errors: 7

## Errors

### Critical (core functionality)

- [ ] `ComIndexMulti()` with multiple indices: crashes with `'list' object cannot be coerced to type 'double'` when more than one index expression is passed -- (tested with `index = c("mean(...)", "sd(...)")`) -- The internal aggregation loop attempts to coerce a list of matrices to a numeric vector; reviewed under issue tracking as a pre-existing bug. Workaround: call `ComIndexMulti()` separately for each index.

### Plotting

- [ ] `plot.ComIndexMulti()`: errors with `'from' must be a finite number` when `ComIndexMulti` results contain all-`NaN` values (e.g. when `by.factor` groups produce empty cells) -- (tested with audit synthetic data where `by.factor == ind.plot`) -- Add a guard that checks for all-`NaN`/all-`NA` SES matrices before calling `plot.listofindex()`, or emit an informative error.

- [ ] `plot.traitflex()` with `use.percentage = FALSE`: crashes with `argument "height" is missing, with no default` -- (tested with `plot.traitflex(res, use.percentage = FALSE)`) -- The `use.percentage` code path calls `barplot()` without supplying the `height` argument. Fix the internal call to pass the correct data matrix.

- [ ] `plotSpVar()` with traits containing `NA`: crashes with `attempt to set 'colnames' on an object with less than two dimensions` -- (tested with a trait matrix containing `NA` values) -- The function drops `NA` rows before checking that the resulting object has at least two dimensions; add a dimension check or `drop = FALSE` guard after `NA` removal.

- [ ] `plotDistri()`: errors with `'legend' is of length 0` -- (tested with standard trait/plot vectors) -- The legend label vector is empty when plot groups have no names; ensure factor levels are used as legend labels with a fallback.

- [ ] `plotSESvar()`: errors with `argument is of length zero` -- (tested with a valid `Tstats` result) -- An internal index or name lookup returns a zero-length vector when trait names are missing or NULL; add a `!is.null` / `length > 0` guard before indexing.

### Input validation

- [ ] `AbToInd()`: crashes with `number of rows of traits and com need to be equal` when the number of species in `traits` does not match the number of rows in `com` -- (tested with mismatched audit data) -- The error message is informative, but could be improved to report actual dimensions. This is expected behaviour for invalid input; no code change required, but the existing stop message could name the arguments.

## Recommendations

- [ ] `ComIndexMulti()`: the multi-index crash (see above) should be treated as a blocking bug before the next CRAN release.
- [ ] Plotting functions (`plotSpVar`, `plotDistri`, `plotSESvar`, `plot.traitflex`): none of these are tested; adding at least one `expect_no_error` test per function would catch regressions early.
- [ ] `ses()` / `Pval()`: the square-obs warning message references "first and second dimension" of the null model, which is now incorrect for the common case where permutations are stored in the first dimension. The warning text should be updated.

## Skipped functions

None — all exported functions were exercised.

## Fixed during this audit cycle

The following bugs were identified during audit preparation and fixed in `R/allfunctions_cati.R`:

- `ses()` / `ses.Tstats()` / `ses.listofindex()` / `plot.Tstats()` / `plot.listofindex()`: "non-conformable arrays" when number of communities equals number of traits (square observed matrix). The null model's permutation dimension is now detected automatically.
- `ses()` / `Pval()` / `plot.ComIndexMulti()`: R >= 4.5 error "the condition has length > 1" caused by `if (class(x) == "list")` where `class()` returns a two-element vector for matrices. Fixed with `is.list()`.
- `sum_Tstats()` type `"percent"` and `"all"`: "NAs are not allowed in subscripted assignments" caused by NA values in logical subscripts. Fixed with `which()` wrappers in both code blocks (the function contained duplicate assignment blocks, both needed patching).
