# Plot result of observed indices values against null distribution

Function to plot result of observed indices values against null
distribution.

## Usage

``` r
plotRandtest(x, alternative = "two-sided", ...)
```

## Arguments

- x:

  An object of class listofindex, ComIndex, ComIndexMulti or Tstats.

- alternative:

  Indicates the alternative hypothesis and must be one of "two.sided",
  "greater" or "less". You can specify just the initial letter.
  "greater" corresponds to positive association, "less" to negative
  association.

- ...:

  Any additional arguments are passed to the plot function creating the
  core of the plot and can be used to adjust the look of resulting
  graph.

## Value

None; used for the side-effect of producing a plot.

## Author

Adrien Taudiere

## See also

[`ComIndex`](https://adrientaudiere.github.io/cati/reference/ComIndex.md);
[`ComIndexMulti`](https://adrientaudiere.github.io/cati/reference/ComIndexMulti.md);
[`Tstats`](https://adrientaudiere.github.io/cati/reference/Tstats.md);
[`as.listofindex`](https://adrientaudiere.github.io/cati/reference/as.listofindex.md);
[`plot.listofindex`](https://adrientaudiere.github.io/cati/reference/plot.listofindex.md)

## Examples

``` r
  data(finch.ind)
  # \donttest{
    res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch,
    sp = sp.finch, nperm = 9, printprogress = FALSE)
#> Warning: This function exclude 1137 Na values
  plotRandtest(res.finch)




















































































  # }
```
