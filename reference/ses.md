# Standardized effect size and confidence interval for a matrix of statistics

calculation standardized effect size and confidence interval for a
matrix of statistics and the related null model expressed as a list or
as an array. Internal function use by other functions of the package.
You can transpose the observed matrix to represent either the SES by
traits or by plots.

## Usage

``` r
ses(obs = NULL, nullmodel = NULL, val.quant = c(0.025, 0.975))
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

## Details

Warning: to detect automatically the correspondence between dimension of
observed matrix and null model list or array, observed matrix needs to
have different numbers of rows and columns. In the case of same row and
column number, please verify manually the correspondance beatween the
rows of the observed matrix and the null model array.

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
#>               WingL      BeakH     UBeakL     N.UBkL
#> DMaj      -2.827416  -1.336514  -2.074583  -6.751185
#> EspHd    -10.068256 -11.006279 -13.611751 -12.631524
#> FlorChrl -21.789018 -25.929804 -20.774729 -32.547878
#> GnovTwr  -59.454754 -56.663400 -39.220716 -58.351808
#> MrchBndl -38.246567 -15.386426 -15.512534 -50.037458
#> SCruInde -15.142589 -22.065564 -20.213692 -52.348914
#> 
#> $ses.inf
#>          WingL      BeakH     UBeakL    N.UBkL
#> [1,] -1.007129 -1.2381827 -1.0524118 -1.563383
#> [2,] -1.393112 -0.9323729 -1.4059753 -1.458427
#> [3,] -1.612186 -1.2769546 -1.1530722 -1.911265
#> [4,] -1.180486 -1.9737183 -1.3388203 -1.215090
#> [5,] -1.312835 -1.4099973 -1.9800838 -1.116999
#> [6,] -1.088765 -1.3031209 -0.8208645 -1.215098
#> 
#> $ses.sup
#>         WingL     BeakH   UBeakL   N.UBkL
#> [1,] 1.657229 1.6205856 1.896485 1.420843
#> [2,] 1.216337 2.0328700 1.503473 1.182347
#> [3,] 1.486941 1.5023207 1.551035 1.212148
#> [4,] 1.919877 0.8223272 1.495129 1.531954
#> [5,] 1.202252 1.2602863 1.068857 1.524495
#> [6,] 1.894855 1.3727849 2.012200 1.734555
#> 
#> attr(,"class")
#> [1] "ses"
  # }
```
