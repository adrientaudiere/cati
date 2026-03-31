# Plot SES values against a variable

Plot standardized effect size values against a variable

## Usage

``` r
plotSESvar(index.list, variable = NULL, ylab = "variable", 
  color.traits = NULL, val.quant = c(0.025, 0.975), resume = FALSE, 
  multipanel = TRUE)
```

## Arguments

- index.list:

  A list of index and the associate null models in the forme: list(
  index_1 = index_1_observed, index_1_nm = null.model.index_1 ,index_2 =
  index_2_observed, index_2_nm = null.model.index_2, ...).

- variable:

  The variable against standardized effect sizes are plotted.

- ylab:

  Label for the variable.

- color.traits:

  A vector of colors corresponding to traits.

- val.quant:

  Numeric vectors of length 2, giving the quantile to calculation
  confidence interval. By default val.quant = c(0.025,0.975) for a
  bilateral test with alpha = 5%.

- resume:

  Logical value; resume = FALSE by default; Simplify the plot by
  plotting the mean and standard error for index value of multiple
  traits

- multipanel:

  Logical value. If TRUE divides the device to shown several traits
  graphics in the same device.

## Value

None; used for the side-effect of producing a plot.

## Author

Adrien Taudiere

## See also

[`plot.listofindex`](https://adrientaudiere.github.io/cati/reference/plot.listofindex.md);
[`ses`](https://adrientaudiere.github.io/cati/reference/ses.md)

## Examples

``` r
  data(finch.ind)
  res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, sp = sp.finch, 
  nperm = 9)
#> Warning: This function exclude 1137 Na values
#> [1] "creating null models"
#> [1] "8.33 %"
#> [1] "16.67 %"
#> [1] "25 %"
#> [1] "33.33 %"
#> [1] "41.63 %"
#> [1] "49.97 %"
#> [1] "58.3 %"
#> [1] "66.63 %"
#> [1] "74.93 %"
#> [1] "83.27 %"
#> [1] "91.6 %"
#> [1] "99.93 %"
#> [1] "calculation of Tstats using null models"
#> [1] "8.33 %"
#> [1] "16.67 %"
#> [1] "25 %"
#> [1] "33.33 %"
#> [1] "41.63 %"
#> [1] "49.97 %"
#> [1] "58.3 %"
#> [1] "66.63 %"
#> [1] "74.93 %"
#> [1] "83.27 %"
#> [1] "91.6 %"
#> [1] "99.93 %"
# \donttest{
  par(mfrow = c(2,2))
  species.richness <- table(ind.plot.finch)
  plotSESvar(as.listofindex(list(res.finch)), species.richness, 
  multipanel = FALSE)


  #Same plot with resume = TRUE.
  
  par(mfrow = c(2,2))
  plotSESvar(as.listofindex(list(res.finch)), species.richness, 
  resume = TRUE, multipanel = FALSE)

  par(mfrow = c(1,1))
  # }
```
