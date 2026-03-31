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
#>               WingL      BeakH    UBeakL     N.UBkL
#> DMaj      -1.870856  -1.213164  -1.54033  -2.822821
#> EspHd    -14.014723 -12.308099 -25.06424 -19.768305
#> FlorChrl -23.447999 -18.319063 -24.39297 -36.404323
#> GnovTwr  -65.842469 -44.462579 -39.57467 -65.106153
#> MrchBndl -36.910612 -13.886818 -20.81403 -31.402890
#> SCruInde -13.623919 -26.365959 -35.17858 -54.975673
#> 
#> $ses.inf
#>          WingL     BeakH    UBeakL    N.UBkL
#> [1,] -1.356619 -1.016641 -1.263042 -1.195067
#> [2,] -1.690393 -1.377036 -1.087942 -1.284178
#> [3,] -1.816127 -1.458809 -1.195481 -1.279466
#> [4,] -0.999796 -1.263704 -1.657674 -2.005284
#> [5,] -1.564175 -1.723526 -1.764895 -1.917702
#> [6,] -1.206963 -1.436822 -1.698174 -1.282041
#> 
#> $ses.sup
#>          WingL    BeakH   UBeakL    N.UBkL
#> [1,] 1.4911218 1.885793 1.763568 1.8444413
#> [2,] 1.0297603 1.611896 1.751019 1.1735537
#> [3,] 0.9957055 1.171949 1.308663 1.6089367
#> [4,] 1.6437235 1.646858 1.255040 0.9807688
#> [5,] 1.3498685 1.254512 1.365385 1.1382232
#> [6,] 1.5765458 1.222468 1.453370 1.2769790
#> 
#> attr(,"class")
#> [1] "ses"
#> 
#> $index_1_2
#> $ses
#>                WingL      BeakH     UBeakL      N.UBkL
#> DMaj      -4.0074226  -6.134675  -4.413250  -7.9497810
#> EspHd      6.0805779   5.171553   2.323167  18.3888596
#> FlorChrl -16.7991154 -29.154401 -22.934417 -24.1316369
#> GnovTwr    9.9665127  17.251637  11.525304   5.7818044
#> MrchBndl   6.3069224   8.527765  -1.663504  -0.2460368
#> SCruInde  -0.3132923  -3.817210  -1.691014  -0.6774189
#> 
#> $ses.inf
#>          WingL     BeakH    UBeakL    N.UBkL
#> [1,] -1.063288 -1.389213 -1.539777 -1.339774
#> [2,] -1.665485 -1.912845 -1.057534 -1.478740
#> [3,] -1.556002 -1.715757 -1.401600 -1.408452
#> [4,] -1.748082 -1.754989 -1.811362 -1.846781
#> [5,] -1.186025 -1.345218 -1.387292 -1.392932
#> [6,] -1.722516 -1.079109 -1.503882 -1.429248
#> 
#> $ses.sup
#>         WingL     BeakH   UBeakL   N.UBkL
#> [1,] 1.324373 1.0713318 1.322398 1.829593
#> [2,] 1.086929 1.0693983 1.791520 1.657671
#> [3,] 1.244269 1.0539682 1.032648 1.539066
#> [4,] 1.400800 0.9789668 1.206622 1.069753
#> [5,] 1.491693 1.2799229 1.317195 1.198347
#> [6,] 1.410351 1.6776525 1.278493 1.441474
#> 
#> attr(,"class")
#> [1] "ses"
#> 
#> $index_1_3
#> $ses
#>               WingL      BeakH     UBeakL     N.UBkL
#> DMaj     -1.0628125         NA  0.2878515 -0.8039819
#> EspHd     1.6758376  1.2516334  0.2739161  1.2849256
#> FlorChrl  1.6463877  0.1745596  0.5552939 -0.8119549
#> GnovTwr   0.4536874  4.2327413  1.8113528  0.6145748
#> MrchBndl  0.7364839  0.1511099 -0.4557055  1.3117143
#> SCruInde  0.2410911 -0.6317473  0.2209433  3.6305713
#> 
#> $ses.inf
#>           WingL     BeakH     UBeakL     N.UBkL
#> [1,] -1.0706749        NA -0.5022323 -0.8754424
#> [2,] -1.0086646 -1.073249 -1.0153685 -0.6450516
#> [3,] -1.1576953 -1.534509 -0.7969440 -1.5118315
#> [4,] -0.9888005 -1.662890 -0.9693731 -1.0776340
#> [5,] -0.9839601 -1.173680 -1.0417689 -0.8615177
#> [6,] -1.2388702 -1.515353 -1.0754667 -1.3713547
#> 
#> $ses.sup
#>         WingL    BeakH   UBeakL   N.UBkL
#> [1,] 1.703741       NA 1.350328 1.546335
#> [2,] 1.655742 1.577536 1.713173 2.090218
#> [3,] 1.606560 1.493578 2.023353 1.144518
#> [4,] 1.707577 1.104840 1.713861 1.551609
#> [5,] 1.742872 1.420183 1.724377 1.818013
#> [6,] 1.694147 1.026320 1.481175 1.520574
#> 
#> attr(,"class")
#> [1] "ses"
#> 
#> attr(,"class")
#> [1] "ses.list"
# }
```
