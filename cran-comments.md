## Submission

This is a bug-fix release (0.99.6 -> 0.99.7). It fixes the bugs reported in an
external audit of 0.99.6, several of which gave wrong numerical results in
`Tstats()`, `ComIndex()`, `ComIndexMulti()`, `ses()` and `Pval()`. NEWS.md lists
every change, including a breaking-changes section (`Pval()` tail names now
match `sum_Tstats()`, `Pval()` returns every trait, `Tstats()` rows follow
`levels(ind.plot)`).

## R CMD check results

Local check with `--as-cran`, `--run-donttest` and remote incoming checks:

0 errors | 0 warnings | 0 notes

## Notes on urlchecker

* `urlchecker::url_check()` timed out on `https://CRAN.R-project.org/package=cati`
  (README.md badges) from the local network. The URL is valid.

## Downstream dependencies

There are no downstream dependencies for this package.
