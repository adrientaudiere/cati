# Ratio of the shortest distance to the longest distance in a minimum spanning tree

Ratio of the shortest distance to the longest distance in a minimum
spanning tree.

## Usage

``` r
MinMaxMST(traits, gower.dist = TRUE, scale.tr = TRUE, method.dist = "euclidian")
```

## Arguments

- traits:

  Traits matrix (traits in column)

- gower.dist:

  Calculate Gower distance using the function `daisy` from package
  cluster.

- scale.tr:

  Does traits need to be scale before multi-traits metric calculation?
  Only use when gower.dist = FALSE. Default is yes.

- method.dist:

  Method to calculate the distance in case of multi-traits metric
  (function dist). Only use when gower.dist = FALSE. Default is
  euclidian.

## Value

The value of the ratio of the shortest distance to the longest distance
in a minimum spanning tree.

## References

Stubbs, WJ., and Wilson, JB. 2004. Evidence for limiting similarity in a
sand dune community. Journal of Ecology 92: 557-567. Aiba, M.,
Katabuchi, M., Takafumi, H., Matsuzaki, S.S., Sasaki, T. & Hiura, T.
2013. Robustness of trait distribution metrics for community assembly
studies under the uncertainties of assembly processes. Ecology, 94,
2873-2885.

## Author

Aiba et al., 2013 modified by Adrien Taudiere

## Examples

``` r
# \donttest{

  data(finch.ind)
  
  MinMaxMST(traits.finch[1:10,])
#> [1] 0.1059942
  MinMaxMST(traits.finch[1:10,], gower.dist = FALSE)
#> [1] 0.1346956
  MinMaxMST(traits.finch[1:10,], gower.dist = FALSE, scale.tr = FALSE)
#> [1] 0.1259882
# }
```
