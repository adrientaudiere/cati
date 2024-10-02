# cati: Communities Assembly by Traits: Individuals and beyond

[![GPLv3](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl.html)
[![Travis-CI Build Status](https://travis-ci.org/adrientaudiere/cati.svg?branch=master)](https://travis-ci.org/adrientaudiere/cati)
[![CRAN by month](http://cranlogs.r-pkg.org/badges/cati?color=red)](http://cran.rstudio.com/web/packages/cati/index.html)
[![CRAN Total](http://cranlogs.r-pkg.org/badges/grand-total/cati?color=yellowgreen)](http://cran.rstudio.com/web/packages/cati/index.html)
[![DOI](https://zenodo.org/badge/19670/adrientaudiere/cati.svg)](https://zenodo.org/badge/latestdoi/19670/adrientaudiere/cati)
[![Research software impact](http://depsy.org/api/package/cran/cati/badge.svg)](http://depsy.org/package/r/cati)

## A R package to detect communities assembly processes by functionnal traits

The package is described in *Ecography*: **[cati: an R package using functional traits to detect and quantify
multi-level community assembly processes](http://onlinelibrary.wiley.com/doi/10.1111/ecog.01433/pdf)**.


### Install the latest versions of all dependencies from CRAN:

```r
install.packages(c("devtools", "e1071", "mice", "rasterVis", "hypervolume", "FD", "geometry", "vegan", "nlme", "ade4", "ape"))
```
### Install cati's current development version from Github:

```r
devtools::install_github("adrientaudiere/cati")
```

### Install cati's cran version:

```r
install.packages("cati")
```

### Attach the package and you are ready to start:
```r
library(cati)
```


There's a **[tutorial](https://github.com/adrientaudiere/cati/blob/Package-cati/Documentation/vignette_Darwin_finches/vignette.pdf)** which illustrate the cati package using Darwin finches data.


## How to cite?

"Taudiere, A. and Violle, C. (2016), cati: an R package using functional traits to detect and quantify multi-level community assembly processes. Ecography. doi: 10.1111/ecog.01433"
