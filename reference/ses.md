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
#>              WingL     BeakH     UBeakL     N.UBkL
#> DMaj      -2.18895  -2.21453  -1.822727  -2.596395
#> EspHd    -16.43843  -7.82619 -27.944877 -21.364184
#> FlorChrl -24.86371 -20.39617 -17.498187 -27.563940
#> GnovTwr  -36.95583 -41.72193 -27.109372 -89.621058
#> MrchBndl -24.19264 -17.65312 -23.215739 -23.315833
#> SCruInde -19.32696 -16.49325 -26.120577 -32.599714
#> 
#> $ses.inf
#>          WingL      BeakH    UBeakL    N.UBkL
#> [1,] -1.049982 -1.1063311 -1.092537 -1.382870
#> [2,] -1.319003 -0.8447874 -1.379337 -1.565618
#> [3,] -1.289835 -1.5857315 -1.479637 -1.468505
#> [4,] -1.326733 -1.2322373 -2.041865 -1.506082
#> [5,] -1.823913 -1.4779954 -2.060887 -1.715350
#> [6,] -1.645602 -1.3731655 -1.339902 -1.259153
#> 
#> $ses.sup
#>         WingL    BeakH    UBeakL    N.UBkL
#> [1,] 1.361774 1.685808 1.5376225 1.4212593
#> [2,] 1.515642 1.826772 1.1954235 1.5130455
#> [3,] 1.702778 1.683091 1.4492612 1.5791214
#> [4,] 1.541997 1.788785 0.8760959 1.0518123
#> [5,] 1.274221 1.527393 0.8046066 0.9911008
#> [6,] 1.410789 1.728268 1.4281265 1.5087821
#> 
#> attr(,"class")
#> [1] "ses"
  # }
```
