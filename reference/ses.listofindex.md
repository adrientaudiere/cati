# Standardized effect size for a list of index.

Standardized effect size and confidence interval for a list of index.

## Usage

``` r
ses.listofindex(index.list = NULL, val.quant = c(0.025, 0.975))
```

## Arguments

- index.list:

  A list of index obtain using the function as.listofindex.

- val.quant:

  Numeric vectors of length 2, giving the quantile to calculation
  confidence interval. By default val.quant = c(0.025,0.975) for a
  bilateral test with alpha = 5%.

## Value

A list which each component correspond to the result of the ses function
for an index. Further, each component is a list of three components:

- \$ses :

  Observed value of standardized effect size.

- \$ses.inf :

  Lower limit of the confidence interval.

- \$ses.sup :

  Upper limit of the confidence interval.

## Author

Adrien Taudiere

## See also

[`as.listofindex`](https://adrientaudiere.github.io/cati/reference/as.listofindex.md);
[`ses`](https://adrientaudiere.github.io/cati/reference/ses.md)

## Examples

``` r
# \donttest{
data(finch.ind)
res.finch <- Tstats(traits.finch, ind.plot = ind.plot.finch,
  sp = sp.finch, nperm = 9, print = FALSE)
#> Warning: This function exclude 1137 Na values
index.list <- as.listofindex(list(res.finch))
ses.listofindex(index.list)
#> $index_1_1
#> $ses
#>               WingL       BeakH     UBeakL     N.UBkL
#> DMaj      -3.752264  -0.8198705  -2.577225  -1.505947
#> EspHd    -30.342358 -15.8767416 -14.456797 -34.025467
#> FlorChrl -18.291903 -21.7330573 -28.264245 -27.465682
#> GnovTwr  -60.679064 -36.0085767 -93.751560 -84.400695
#> MrchBndl -32.487132 -19.5580333 -30.629497 -38.005418
#> SCruInde -23.177564 -20.1618509 -44.861638 -25.882607
#> 
#> $ses.inf
#>          WingL      BeakH    UBeakL     N.UBkL
#> [1,] -1.591954 -0.7888253 -1.108731 -1.1180194
#> [2,] -1.117031 -1.8915128 -1.116845 -1.6492181
#> [3,] -1.264779 -1.1952185 -1.450756 -0.9954512
#> [4,] -1.612978 -1.3570750 -1.401471 -2.0590084
#> [5,] -1.527233 -1.7332233 -1.041122 -1.1794593
#> [6,] -1.180212 -0.7554718 -1.085056 -1.4773447
#> 
#> $ses.sup
#>         WingL    BeakH   UBeakL    N.UBkL
#> [1,] 1.294122 1.936475 1.689619 1.5394236
#> [2,] 1.463884 1.194471 1.305519 1.1711942
#> [3,] 1.678640 1.323098 1.707068 1.9638947
#> [4,] 1.225489 1.448174 1.137625 0.7812745
#> [5,] 1.581308 1.114030 1.894577 1.6998856
#> [6,] 1.500286 1.904737 1.569963 1.2679092
#> 
#> attr(,"class")
#> [1] "ses"
#> 
#> $index_1_2
#> $ses
#>               WingL      BeakH      UBeakL      N.UBkL
#> DMaj      -3.736519  -3.075777  -3.1165419  -4.7761695
#> EspHd     14.929840   5.691378   1.5230704   7.6521738
#> FlorChrl -31.563685 -15.208449 -12.9882451 -16.2719218
#> GnovTwr   13.723690  11.332534  15.4103303  12.0869920
#> MrchBndl   5.214631  10.098249  -1.9488387  -0.7946450
#> SCruInde  -2.532138  -5.997187  -0.8944711  -0.3054629
#> 
#> $ses.inf
#>           WingL     BeakH    UBeakL    N.UBkL
#> [1,] -1.6040604 -1.559732 -1.694023 -1.574192
#> [2,] -1.2386301 -1.216565 -1.756959 -1.685485
#> [3,] -1.3636885 -1.050502 -1.484664 -1.539251
#> [4,] -1.3858760 -1.531473 -1.193741 -1.392077
#> [5,] -0.9032113 -1.461069 -1.217765 -1.020906
#> [6,] -1.3671906 -1.371484 -1.698442 -1.674502
#> 
#> $ses.sup
#>         WingL    BeakH   UBeakL   N.UBkL
#> [1,] 1.428656 1.638703 1.252801 1.648485
#> [2,] 1.673489 1.442731 1.143902 1.424790
#> [3,] 1.578594 1.629177 1.376524 1.511298
#> [4,] 1.399090 1.251809 1.494159 1.340665
#> [5,] 1.875025 1.472683 1.722506 1.554509
#> [6,] 1.578073 1.039978 1.309521 1.225101
#> 
#> attr(,"class")
#> [1] "ses"
#> 
#> $index_1_3
#> $ses
#>                WingL      BeakH     UBeakL     N.UBkL
#> DMaj     -0.89625014 -0.9997452 -0.4457254 -0.8205937
#> EspHd     9.20814581  2.4123278  3.3402890 13.9525127
#> FlorChrl -1.12295901 -1.4614574 -1.1180587 -1.5563874
#> GnovTwr  -0.06040716  1.2497877  0.4201516  0.9604926
#> MrchBndl 10.87603175  1.4942136  1.3377886  6.4287902
#> SCruInde  7.18170085  6.1492987  2.9937989  7.6138571
#> 
#> $ses.inf
#>           WingL      BeakH     UBeakL     N.UBkL
#> [1,] -1.0191991 -0.9019358 -0.7888810 -1.3045226
#> [2,] -1.5410105 -0.6816131 -0.7409673 -0.8270074
#> [3,] -1.6334123 -1.5704769 -1.6403229 -1.8104180
#> [4,] -0.8163395 -0.6788379 -0.9981933 -1.1719520
#> [5,] -1.2952494 -0.9887366 -0.6782023 -0.9265183
#> [6,] -1.2240411 -1.1747463 -0.7940310 -1.3515922
#> 
#> $ses.sup
#>         WingL    BeakH   UBeakL    N.UBkL
#> [1,] 1.732765 2.010743 1.971815 1.6617608
#> [2,] 1.462021 2.090148 1.989423 1.7340187
#> [3,] 1.140418 1.012266 1.102982 0.8917151
#> [4,] 1.891222 2.056071 1.512372 1.6038961
#> [5,] 1.797869 1.904356 1.944364 1.5745008
#> [6,] 1.359412 1.819281 1.982604 1.6113863
#> 
#> attr(,"class")
#> [1] "ses"
#> 
#> attr(,"class")
#> [1] "ses.list"
# }
```
