# cati: Communities Assembly by Traits: Individuals and beyond

[![Travis-CI Build Status](https://travis-ci.org/adrientaudiere/cati.svg?branch=master)](https://travis-ci.org/adrientaudiere/cati)
[![Downloads from the RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/cati)](http://cran.rstudio.com/package=cati)
[![Coverage Status](https://coveralls.io/repos/github/adrientaudiere/cati/badge.svg?branch=master)](https://coveralls.io/github/adrientaudiere/cati?branch=master)
[![DOI](https://zenodo.org/badge/19670/adrientaudiere/cati.svg)](https://zenodo.org/badge/latestdoi/19670/adrientaudiere/cati)

---
## A R package to detect communities assembly processes by functionnal traits

The package is described in *Ecography*: **[cati: an R package using functional traits to detect and quantify
multi-level community assembly processes](http://onlinelibrary.wiley.com/doi/10.1111/ecog.01433/pdf)**.


### Install the latest versions of all dependencies from CRAN:

```r
install.packages(c("devtools", "e1071", "mice", "rasterVis", "hypervolume", "FD", "geometry", "vegan", "nlme", "ade4", "ape"))
```
### Install cati's current development version from Github:

```r
devtools::install_github("adrientaudiere/cati/pkg/cati")
```

### Install cati's cran version:
```r
install.packages("cati")
```

### Attach the package and you are ready to start:
```r
library(cati)
```

---
There's a **[tutorial](https://github.com/adrientaudiere/cati/blob/Package-cati/Documentation/vignette_Darwin_finches/vignette.pdf)** which illustrate the cati package using Darwin finches data.

---
How to cite?

"Taudiere, A. and Violle, C. (2015), cati: an R package using functional traits to detect and quantify multi-level community assembly processes. Ecography. doi: 10.1111/ecog.01433"
