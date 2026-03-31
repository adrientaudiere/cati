# Coefficient of variation, mean, minimum and standard deviation of the nearest neigbourhood distance.

CVNND : Coefficient of variation of the nearest neigbourhood distance

MNND : Mean of the nearest neigbourhood distance

MinNND : Minimum of the nearest neigbourhood distance

SDNND : Standard deviation of the nearest neigbourhood distance

SDND : Standard deviation of the neigbourhood distance

MND : Mean of the neigbourhood distance

## Usage

``` r
CVNND(traits, div_range =  FALSE, na.rm = FALSE, scale.tr = TRUE,
  method.dist = "euclidian")
  
  MNND(traits, div_range =  FALSE, na.rm = FALSE, scale.tr = TRUE, 
  method.dist = "euclidian")
  
  MinNND(traits, div_range =  FALSE, na.rm = FALSE, scale.tr = TRUE, 
  method.dist = "euclidian")
  
  SDNND(traits, div_range =  FALSE, na.rm = FALSE, scale.tr = TRUE, 
  method.dist = "euclidian")
  
  SDND(trait, div_range = FALSE, na.rm = FALSE)  
  
  MND(trait, div_range = FALSE, na.rm = FALSE)
```

## Arguments

- traits:

  Trait vector (uni-trait metric) or traits matrix (Multi-traits
  metric), traits in column.

- trait:

  Trait vector

- div_range:

  Does metric need to be divided by the range? Default is no.

- na.rm:

  If div_range=TRUE, a logical value indicating whether NA values should
  be stripped before the computation proceeds.

- scale.tr:

  Does traits need to be scale before multi-traits metric calculation?
  Default is yes.

- method.dist:

  Method to calculate the distance in case of multi-traits metric
  (function dist). Default is euclidian.

## Value

One value corresponding to the metric value.

## References

Aiba, M., Katabuchi, M., Takafumi, H., Matsuzaki, S.S., Sasaki, T. &
Hiura, T. 2013. Robustness of trait distribution metrics for community
assembly studies under the uncertainties of assembly processes. Ecology,
94, 2873-2885. Jung, Vincent, Cyrille Violle, Cedric Mondy, Lucien
Hoffmann, et Serge Muller. 2010. Intraspecific variability and
trait-based community assembly: Intraspecific variability and community
assembly. Journal of Ecology 98 (5): 1134-1140.

## Author

Adrien Taudiere

## Examples

``` r
data(finch.ind)

CVNND(traits.finch[,1], na.rm = TRUE)
#> [1] 49.62862
CVNND(traits.finch[,1], div_range =  TRUE, na.rm = TRUE)
#> [1] 10.65975
CVNND(traits.finch, na.rm = TRUE)
#> [1] 0.6326296
CVNND(traits.finch, scale.tr = FALSE, na.rm = TRUE)
#> [1] 0.7362636
SDND(traits.finch[,1], na.rm = TRUE)
#> [1] 0.1248847
```
