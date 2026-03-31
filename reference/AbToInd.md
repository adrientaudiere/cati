# Internal function. Transform abundance data matrix into individual like matrix.

Transform abundance data matrix into individual like matrix to allows
the use of ComIndex and ComIndexMulti on populationnal or specific
traits values.

## Usage

``` r
AbToInd(traits, com, type.sp.val = "count")
```

## Arguments

- traits:

  Individual Matrix of traits with traits in columns. "traits" matrix
  must have row names (e.g. species or populationnal names).

- com:

  Community data matrix with species in rows and sites in column.

- type.sp.val:

  Either "count" or "abundance". Use abundance when all values in the
  com matrix are not superior to one. Using abundance is EXPERIMENTAL.
  This function round abundance to fit count data.

## Details

Internal function

## Value

A list of objects:

- \$traits:

  Individual traits matrix

- \$sp:

  Vector of species attributes

- \$ind.plot:

  Vector of sites attributes

## Author

Adrien Taudiere
