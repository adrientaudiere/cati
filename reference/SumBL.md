# Sum of branch length of a classification dendrogram (Petchey and Gaston, 2002)

Sum of branch length of a classification dendrogram (Petchey and Gaston,
2002)

## Usage

``` r
SumBL(traits, gower.dist = TRUE, method.hclust = "average", 
  scale.tr = TRUE, method.dist = "euclidian")
```

## Arguments

- traits:

  Traits matrix (traits in column)

- gower.dist:

  Calculate Gower distance using the function `daisy` from package
  cluster.

- method.hclust:

  Define the method for the hclust function (default is "average" i.e.
  UPGMA).

- scale.tr:

  Does traits need to be scale before multi-traits metric calculation?
  Only use when gower.dist = FALSE. Default is yes.

- method.dist:

  Method to calculate the distance in case of multi-traits metric
  (function dist). Only use when gower.dist = FALSE. Default is
  euclidian.

## Value

The value of the sum of branch length from a classification dendrogram
of traits.

## References

Petchey, OL., and Gaston, KJ. 2002. Functional diversity (FD), species
richness and community composition. Ecology Letters 5:402-411

## Author

Adrien Taudiere

## Examples

``` r
# \donttest{

data(finch.ind)
SumBL(traits.finch)
#> Warning: All individuals with one NA are excluded
#> [1] 28.98831
SumBL(traits.finch, gower.dist = FALSE)
#> Warning: All individuals with one NA are excluded
#> [1] 345.3096

# }
```
