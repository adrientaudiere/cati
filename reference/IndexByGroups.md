# Apply metrics to groups.

Transforme a list of metrics to apply them to groups, typically to
populations.

## Usage

``` r
IndexByGroups(metrics, groups)
```

## Arguments

- metrics:

  A vector of metrics like the argument "index" of function ComIndex

- groups:

  Name of the factor to apply the metrics to groups in the form "pop",
  e.g. population

## Value

A vector of transformed metrics

## Author

Adrien Taudiere

## Examples

``` r
IndexByGroups(c("mean(x)", "sd(x)"), "pop")
#> [1] "tapply(x, pop, function(x) mean(x))" "tapply(x, pop, function(x) sd(x))"  
```
