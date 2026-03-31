# Functional richness, evenness and divergence following Villeger et al. 2008

Compute the 3 functional diversity indices (multi-traits) presented in
Villeger et al. 2008 (Ecology 89 2290-2301): Functional richness (FRic),
Functional evenness (FEve), Functional divergence (FDiv)

## Usage

``` r
Fred(traits, ind.plot)
```

## Arguments

- traits:

  Individual Matrix of traits with traits in columns. NA are not allowed
  .

- ind.plot:

  Factor defining the name of the plot in which the individual is.

## Value

list of 4 vectors with values of indices in each sites

- \$nbind:

  number of individuals

- \$FRic:

  functional richness index

- \$FEve:

  functional evenness index

- \$FDiv:

  functional divergence index

## Details

For each trait, values are standardized (mean=0 and standard
deviation=1) For FRic computation, number of individuals must be higher
than number of traits

## Author

Sebastien Villeger sligthy modified by Adrien Taudiere

## See also

[`ComIndexMulti`](https://adrientaudiere.github.io/cati/reference/ComIndexMulti.md)
[`ComIndex`](https://adrientaudiere.github.io/cati/reference/ComIndex.md)

## Examples

``` r
data(finch.ind)
# \donttest{
#For most multivariate functions we need to replace (or exclude) NA values.
#For this example, we use the package mice to complete the data.
if (requireNamespace("mice", quietly = TRUE)) {
  mice_imp <- mice::mice(traits.finch, printFlag = FALSE)
  traits.finch.mice <- mice::complete(mice_imp)
  fred <- Fred(traits.finch.mice, ind.plot.finch)
}
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
# }
```
