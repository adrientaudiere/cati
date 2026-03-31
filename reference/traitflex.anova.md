# Variance decomposition for a given trait used in decompCTRE

This function decomposes variation of trait values within a community
into three sources: (i) the intraspecific trait variability, (ii) the
variability due to species turnover and (iii) their covariation is also
separated. This decomposition is computed for the whole variation in the
trait values and, The formula specified, across the contribution of
various explanatory variables considered in the model. S3 method plot
summarizes graphically the decomposition of trait variation, obtained
with the traitflex.anova function. Print is an other S3 method for
object of class traitflex.

## Usage

``` r
traitflex.anova(formula, specif.avg, const.avg, ...)
  # S3 method for class 'traitflex'
plot(x, plot.total =  FALSE, use.percentage = TRUE, 
  plot.covar =  FALSE, cumul =  FALSE, 
  legend.pos = if (plot.total) "topleft" else "topright", 
  plot.res = TRUE, ...)
  
  # S3 method for class 'traitflex'
print(x, ...)
```

## Arguments

- formula:

  The formula parameter must be a one-sided formula, i.e. starting with
  a tilde character. The response variable is specified by the next two
  arguments, specif.avg and const.avg.

- specif.avg:

  Vector with community trait composition values for a single trait. It
  is calculated from trait values specific to each community (i.e. trait
  values for individual species are 'specific' to each plot, or habitat,
  where the species is found)

- const.avg:

  Vector with community trait composition values for a single trait. It
  is calculated from average (fixed) trait values of individual species
  (i.e. fixed trait value for individual species used for all habitats
  where the species is found)

- x:

  An object of class traitflex.

- plot.total:

  Logical value; if TRUE plot not only the individual components of
  variation, but also the total variation. This is useful particularly
  when the decomposition was done with non-trivial formula (i.e. with
  explanatory variables)

- use.percentage:

  Logical value; if TRUE the individual plotted sources of trait
  variation are shown as percentages of the total variation, on 0-100
  scale.

- plot.covar:

  Logical value; if TRUE the covariance between within-species trait
  variability and the variability due to species composition turnover is
  plotted as yet another category within the stacked bars. The
  plot.covar argument is entirely ignored when plotting traitflex object
  fitted with a formula without any predictor variables.

- cumul:

  Logical value; if TRUE values are shown in a cumulative way.

- legend.pos:

  This argument allows you to specify the position of graph legend. Thus
  argument is entirely ignored when plotting traitflex object created
  with a formula without predictors

- plot.res:

  Logical value; if resume = FALSE plot is not shown but the table of
  values used to print the plot is return.

- ...:

  Optional additional arguments.

## Details

The formula parameter must be a one-sided formula, i.e. starting with a
tilde character. The response variable is specified by the next two
arguments, specif.avg and const.avg.

## Value

An object of class traitflex. There are print and plot methods available
for it. The object contains decomposition of sum of squares into
intraspecific variation component, compositional variation component,
their covariation and total in a SumSq element. This is a data frame
with multiple rows if predictors were specified in formula argument. The
RelSumSq element contains the same table relativized to unit row totals.
Finally, the anova.turnover, anova.total, and anova.diff elements
contain the three aov objects used to decompose the variation.

## References

Leps, Jan, Francesco de Bello, Petr Smilauer and Jiri Dolezal. 2011.
Community trait response to environment: disentangling species turnover
vs intraspecific trait variability effects. Ecography 34 (5): 856-863.

## Author

Jan Leps et al., 2011 modified by Adrien Taudiere

## See also

`print.traitflex`; `plot.traitflex`;
[`decompCTRE`](https://adrientaudiere.github.io/cati/reference/decompCTRE.md)

## Examples

``` r
# \donttest{
data(finch.ind)
res.decomp <- decompCTRE(traits = traits.finch, sp = sp.finch,
  ind.plot = ind.plot.finch, printprogress = FALSE)
#> Warning: All individuals with one NA ( 50 individual value(s)) are excluded for the trait 1 : WingL
#> Warning: All individuals with one NA ( 396 individual value(s)) are excluded for the trait 2 : BeakH
#> Warning: All individuals with one NA ( 668 individual value(s)) are excluded for the trait 3 : UBeakL
#> Warning: All individuals with one NA ( 23 individual value(s)) are excluded for the trait 4 : N.UBkL
tf <- res.decomp[[1]]
print(tf)
#> 
#>  Decomposing trait sum of squares into composition turnover
#>  effect, intraspecific trait variability, and their covariation:
#>       Turnover Intraspec. Covariation  Total
#> Total   4.4935     8.3772     -11.268 1.6022
#> 
#>  Relative contributions:
#>       Turnover Intraspec. Covariation Total
#> Total    2.805      5.228      -7.033     1
plot(tf)

# }
```
