# cati: Communities Assembly by Traits: Individuals and beyond

[![GPLv3](https://www.gnu.org/graphics/gplv3-127x51.png)](https://www.gnu.org/licenses/gpl.html)
[![CRAN by month](https://cranlogs.r-pkg.org/badges/cati?color=red)](https://CRAN.R-project.org/package=cati)
[![CRAN Total](https://cranlogs.r-pkg.org/badges/grand-total/cati?color=yellowgreen)](https://CRAN.R-project.org/package=cati)
[![DOI](https://zenodo.org/badge/19670/adrientaudiere/cati.svg)](https://zenodo.org/badge/latestdoi/19670/adrientaudiere/cati)

## A R package to detect communities assembly processes by functionnal traits

The package is described in *Ecography*: **[cati: an R package using functional traits to detect and quantify multi-level community assembly processes](https://doi.org/10.1111/ecog.01433)**.


### Install the latest versions of all dependencies from CRAN:

```r
install.packages(c("devtools", "e1071", "mice", "rasterVis", "hypervolume", "cluster", "geometry", "vegan", "nlme", "ade4", "ape"))
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
