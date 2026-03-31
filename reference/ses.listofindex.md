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
#>               WingL      BeakH     UBeakL     N.UBkL
#> DMaj      -1.611635  -1.236873  -2.484558  -3.750675
#> EspHd    -14.034912  -6.295014 -21.747111 -15.426916
#> FlorChrl -17.886990 -26.624223 -18.457032 -21.653084
#> GnovTwr  -45.085381 -43.194983 -32.700957 -75.981083
#> MrchBndl -32.655030 -33.197028 -21.985797 -30.593948
#> SCruInde -23.070261 -19.495095 -31.745680 -46.753259
#> 
#> $ses.inf
#>          WingL     BeakH    UBeakL    N.UBkL
#> [1,] -1.195048 -1.141620 -1.315233 -1.783044
#> [2,] -1.238681 -1.049774 -1.141105 -1.554190
#> [3,] -1.558394 -1.801822 -1.355612 -1.865627
#> [4,] -1.700457 -1.302747 -1.928211 -1.338081
#> [5,] -0.786997 -1.370627 -1.730942 -1.839607
#> [6,] -1.478209 -1.275694 -1.207095 -1.554613
#> 
#> $ses.sup
#>         WingL    BeakH    UBeakL    N.UBkL
#> [1,] 1.490762 1.862954 1.5265017 1.2255148
#> [2,] 1.726885 1.789198 1.9312722 1.1129012
#> [3,] 1.665960 1.226424 1.6750862 0.9860102
#> [4,] 1.263595 1.036085 0.9433971 1.3043903
#> [5,] 2.080724 1.643905 0.9085834 1.0439234
#> [6,] 1.297865 1.544098 1.5055888 1.2661495
#> 
#> attr(,"class")
#> [1] "ses"
#> 
#> $index_1_2
#> $ses
#>               WingL      BeakH      UBeakL      N.UBkL
#> DMaj      -4.186735  -3.201492  -4.1116142  -5.9881799
#> EspHd     15.421723   8.245579   3.2828536   7.5180517
#> FlorChrl -15.286903 -14.165286 -16.0902283 -20.1394979
#> GnovTwr   13.242304  14.984858   6.2047368   9.6557613
#> MrchBndl   5.798054   7.298600  -0.9493127   0.1653836
#> SCruInde  -1.110882  -5.270691  -1.1410861   0.1469535
#> 
#> $ses.inf
#>          WingL     BeakH    UBeakL    N.UBkL
#> [1,] -1.190578 -1.437029 -1.345867 -1.225924
#> [2,] -1.242647 -1.177894 -1.230935 -1.332635
#> [3,] -1.791096 -1.603407 -1.022150 -1.039266
#> [4,] -1.337663 -1.334963 -1.528660 -1.279511
#> [5,] -1.459601 -1.362277 -1.372127 -1.734267
#> [6,] -1.246962 -1.621575 -1.411758 -1.422923
#> 
#> $ses.sup
#>         WingL    BeakH   UBeakL   N.UBkL
#> [1,] 1.307677 1.711979 1.174270 1.823600
#> [2,] 1.568678 1.702050 1.780123 1.493693
#> [3,] 1.361643 1.063426 1.504402 1.688133
#> [4,] 1.409762 1.340133 1.385702 1.644915
#> [5,] 1.620050 1.127524 1.633348 1.463988
#> [6,] 1.566848 1.350429 1.513026 1.556173
#> 
#> attr(,"class")
#> [1] "ses"
#> 
#> $index_1_3
#> $ses
#>               WingL      BeakH      UBeakL     N.UBkL
#> DMaj     -0.6947590        Inf -0.51073692 -0.5822476
#> EspHd     1.5056797 -0.1684903  0.23714845  0.3757881
#> FlorChrl  0.3510852  0.9826734 -0.03928972 -0.7216100
#> GnovTwr   1.0607828  2.3334458  0.69723188  0.2962490
#> MrchBndl  0.4707089 -0.5957030 -0.09602140  2.6697564
#> SCruInde -0.1492500 -0.1066701  0.64525200  1.0215812
#> 
#> $ses.inf
#>           WingL      BeakH     UBeakL     N.UBkL
#> [1,] -0.6678087        NaN -0.6717514 -0.9550827
#> [2,] -1.0263043 -0.9395546 -0.5495736 -0.5330648
#> [3,] -0.6822701 -1.6461555 -0.9614212 -1.1798876
#> [4,] -0.4885715 -0.9842290 -0.5931879 -0.7679889
#> [5,] -1.0102708 -1.4479278 -1.0063786 -1.2522417
#> [6,] -1.4975341 -1.2866290 -0.7962683 -0.7708113
#> 
#> $ses.sup
#>         WingL     BeakH    UBeakL   N.UBkL
#> [1,] 2.091366       NaN 0.6717514 1.819760
#> [2,] 1.751572 1.7772762 2.0975206 2.092053
#> [3,] 2.043408 0.9142191 1.3307427 1.600600
#> [4,] 2.111047 1.6198321 2.0777296 1.693805
#> [5,] 1.848137 1.4654305 1.8323221 1.674742
#> [6,] 1.443435 1.3419038 2.0587969 2.016661
#> 
#> attr(,"class")
#> [1] "ses"
#> 
#> attr(,"class")
#> [1] "ses.list"
# }
```
