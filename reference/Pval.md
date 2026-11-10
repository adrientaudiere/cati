# Calcul of p-value for object of class Tstats, ComIndex, ComIndexMulti and listofindex

Calcul of p-value for object of class Tstats, ComIndex, ComIndexMulti
and listofindex. This test equates to finding the quantile in exp in
which obs would be found (under a one-tailed test).

## Usage

``` r
Pval(x, na.rm = TRUE)
```

## Arguments

- x:

  An object of class Tstats, ComIndex, ComIndexMulti or listofindex.
  Objects made by these functions or by
  [`as.listofindex`](https://adrientaudiere.github.io/cati/reference/as.listofindex.md)
  carry the layout of their null values. For a hand-built list whose
  numbers of plots, traits or permutations are equal, the layout is
  guessed from sizes and may be wrong; use
  [`ses`](https://adrientaudiere.github.io/cati/reference/ses.md) with
  `null.dims` in that case.

- na.rm:

  A logical value indicating whether NA values should be stripped before
  the computation proceeds.

## Value

A list of p-value for each metrics, traits and grouping if needed (e.g.
sites). `pval.inf` is small when the observed value is lower than the
null values, `pval.sup` when it is higher. Ties count against
significance; NA null values are dropped; the p-value is NA when the
observed value or all null values are NA.

## Author

Adrien Taudiere

## Examples

``` r

data(finch.ind)
res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
  sp = sp.finch, nperm = 9, print = FALSE)
#> Warning: This function exclude 1137 Na values
  if (FALSE) { # \dontrun{
    Pval(res.finch)
  } # }
```
