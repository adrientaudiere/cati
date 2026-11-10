# Standardized effect size and confidence interval for a matrix of statistics

calculation standardized effect size and confidence interval for a
matrix of statistics and the related null model expressed as a list or
as an array. Internal function use by other functions of the package.
You can transpose the observed matrix to represent either the SES by
traits or by plots; the result is transposed in the same way, except for
a square matrix, whose transposition cannot be detected.

## Usage

``` r
ses(obs = NULL, nullmodel = NULL, val.quant = c(0.025, 0.975), null.dims = NULL)
```

## Arguments

- obs:

  Observed matrix or vector of values.

- nullmodel:

  Either a list or an array of three (two for a vector of observed
  values) dimensions corresponding to the null model permutations.

- val.quant:

  Numeric vectors of length 2, giving the quantile to calculation
  confidence interval. By default val.quant = c(0.025,0.975) for a
  bilateral test with alpha = 5%.

- null.dims:

  Integer vector of length 2 giving the dimensions of the `nullmodel`
  array that match the rows and the columns of `obs`; the remaining
  dimension holds the permutations. For example `c(3, 2)` for a null
  array of permutations x traits x plots and an observed matrix of plots
  x traits. By default (NULL), the layout stored by
  [`as.listofindex`](https://adrientaudiere.github.io/cati/reference/as.listofindex.md)
  is used, or the dimensions are guessed from their sizes.

## Details

Warning: to detect automatically the correspondence between dimension of
observed matrix and null model list or array, observed matrix, null
model and permutations need to have different sizes. When two sizes are
equal, give the layout with `null.dims`. Objects converted with
[`as.listofindex`](https://adrientaudiere.github.io/cati/reference/as.listofindex.md)
(as done by `ses.listofindex`, `Pval` and the plot functions) carry
their layout and need no `null.dims`.

## Value

A list of three components:

- \$ses:

  Observed value of standardized effect size.

- \$ses.inf :

  Lower limit of the confidence interval.

- \$ses.sup :

  Upper limit of the confidence interval.

## Author

Adrien Taudiere

## See also

[`plot.listofindex`](https://adrientaudiere.github.io/cati/reference/plot.listofindex.md);
[`plotSESvar`](https://adrientaudiere.github.io/cati/reference/plotSESvar.md);
[`ses.listofindex`](https://adrientaudiere.github.io/cati/reference/ses.listofindex.md)

## Examples

``` r
  data(finch.ind)
  
  res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch, 
  sp = sp.finch, nperm = 9)
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
    ses(res.finch$Tstats$T_IP.IC, res.finch$Tstats$T_IP.IC_nm)
#> $ses
#>              WingL      BeakH     UBeakL     N.UBkL
#> DMaj      -1.76691  -0.903734  -1.167335  -1.156339
#> EspHd    -15.84097 -17.859382 -13.095555 -13.995191
#> FlorChrl -18.45305 -18.815011 -28.493314 -23.606919
#> GnovTwr  -56.59361 -35.355044 -38.500926 -62.313121
#> MrchBndl -26.19813 -40.636798 -31.019195 -34.914073
#> SCruInde -20.82683 -15.783407 -21.487047 -31.636154
#> 
#> $ses.inf
#>          WingL      BeakH     UBeakL     N.UBkL
#> [1,] -1.324734 -0.9374337 -0.8486808 -0.7224897
#> [2,] -1.452535 -1.2248389 -1.2746699 -1.4012941
#> [3,] -1.103680 -0.9775635 -1.4152969 -0.7871146
#> [4,] -1.148153 -1.4881243 -1.7508124 -1.6390433
#> [5,] -1.618993 -0.9108882 -1.6348139 -1.5728830
#> [6,] -1.285682 -0.9874584 -1.5181695 -1.4771161
#> 
#> $ses.sup
#>         WingL    BeakH   UBeakL   N.UBkL
#> [1,] 1.354442 1.908563 1.768731 1.932252
#> [2,] 1.265432 1.650857 1.330751 1.793144
#> [3,] 1.787032 1.821697 1.802264 2.054830
#> [4,] 1.442650 1.615314 1.203176 1.119197
#> [5,] 1.588414 1.319962 1.480829 1.450005
#> [6,] 1.277500 1.733321 1.173477 1.129879
#> 
#> attr(,"class")
#> [1] "ses"
  # }
```
