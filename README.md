[![DOI](https://zenodo.org/badge/23397/pierreroudier/spemd.svg)](https://zenodo.org/badge/latestdoi/23397/pierreroudier/spemd)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/spemd)](https://cran.r-project.org/package=spemd)

# `spemd`: A bi-dimensional implementation of the Empirical Mode Decomposition (EMD)

## Installation

The package is now available on CRAN:

```
install.packages("spemd")
```

Alternatively, you can install the development version using the `devtools` package:

```
devtools::install_github("pierreroudier/spemd")
```

If you do not have the `devtools` package installed, you can install it simply using:

```
install.packages("devtools")
```

## Quick example

### Gridded data

```
library(spemd)

library(gstat) # to get sample data
library(raster) # for visualisation

# Load sample data
data(ncp.grid, package = 'gstat')
coordinates(ncp.grid) <- ~x+y
gridded(ncp.grid) <- TRUE

plot(raster(ncp.grid))

# Run EMD
res.ncp <- spEMD(ncp.grid, zcol = "depth", thresh.extrema = 0.1, verbose = FALSE)

# Plot results
res <- stack(res.ncp)
plot(res)
```
