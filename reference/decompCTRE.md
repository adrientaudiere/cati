# Variance partitioning for multiple traits

This function decomposes the variation in community trait composition
into three sources: (i) the intraspecific trait variability, (ii) the
variability due to species turnover and (iii) their covariation is also
separated. This decomposition is computed for the whole variation in the
trait values and, The formula specified, across the contribution of
various explanatory variables considered in the model.
Barplot.decompCTRE allow to plot the result of the decomposition.

## Usage

``` r
decompCTRE(traits = NULL, formula = ~1, ind.plot = NULL, sp = NULL, 
  printprogress = TRUE, ...)
  
  # S3 method for class 'decompCTRE'
barplot(height, resume = TRUE, ...)
```

## Arguments

- traits:

  Matrix of traits with traits in column

- height:

  An object of class decompCTRE obtain by the function decompCTRE.

- formula:

  The formula parameter must be a one-sided formula, i.e. starting with
  a tilde (~) character. The response variable is specified by the next
  two arguments, specif.avg and const.avg. By default set to ~1.

- ind.plot:

  Factor defining the name of the plot (site or community) in which the
  individual is.

- sp:

  Factor defining the species which the individual belong to.

- printprogress:

  Logical value; print progress during the calculation or not.

- resume:

  Logical. If resume = FALSE, plot one graphic by traits.

- ...:

  Optional additional arguments

## Value

An object of class "decompCTRE".

## References

Leps, Jan, Francesco de Bello, Petr Smilauer and Jiri Dolezal. 2011.
Community trait response to environment: disentangling species turnover
vs intraspecific trait variability effects. Ecography 34 (5): 856-863.

## Author

Adrien Taudiere Jan Leps

## See also

[`barplot`](https://rdrr.io/r/graphics/barplot.html);
[`traitflex.anova`](https://adrientaudiere.github.io/cati/reference/traitflex.anova.md)

## Examples

``` r
data(finch.ind)
# \donttest{
  res.decomp <- decompCTRE(traits = traits.finch, sp = sp.finch, 
  ind.plot = ind.plot.finch, print = FALSE)
#> Warning: All individuals with one NA ( 50 individual value(s)) are excluded for the trait 1 : WingL
#> Warning: All individuals with one NA ( 396 individual value(s)) are excluded for the trait 2 : BeakH
#> Warning: All individuals with one NA ( 668 individual value(s)) are excluded for the trait 3 : UBeakL
#> Warning: All individuals with one NA ( 23 individual value(s)) are excluded for the trait 4 : N.UBkL

  barplot(res.decomp)


  par(mfrow = c(2,2))
  barplot(res.decomp, resume = FALSE)

  par(mfrow = c(1,1))
  # }
```
