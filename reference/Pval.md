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

- na.rm:

  A logical value indicating whether NA values should be stripped before
  the computation proceeds.

## Value

A list of p-value for each metrics, traits and grouping if needed (e.g.
sites)

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
