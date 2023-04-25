
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Supplementary Materials for: *“Penalized rank-one approximation to functional partial least-squares regression over complex domains”*

> Authors: Harold A. Hernández-Roig (<haroldantonio.hernandez@uc3m.es>);
> M. Carmen Aguilera-Morillo; Eleonora Arnone; Rosa E. Lillo; Laura M.
> Sangalli.

## Overview

- Folder **Animated Betas** contains 3D animations of the estimated
  coefficient functions from *Section 5: “Application to brain
  connectivity data”*.

- Folders **1D_simulations** and **2D_simulations** contains the scripts
  required to reproduce the simulations in *Section 4: “Simulation
  studies”*.

## Simulation studies

The simulations depend on the R package `penR1FPLS`. It can be installed
as follows:

``` r
devtools::install_github("hhroig/penR1FPLS", dependencies = TRUE)
```

### One-dimensional domain

Run:

``` r
source("1D_simulations/main_1D.R")
```

to reproduce the simulations for data defined over a one-dimensional
(1D) domain. The script creates a new folder with results, including
plots.

### Two-dimensional planar domain

Run:

``` r
source("2D_simulations/main_2D.R")
```

to reproduce the simulations for data defined over a two-dimensional
(2D) planar domain. The script creates a new folder with results,
including plots.

### Parallel computing

By default the simulations run in parallel using the `doParallel`
package. For example, in the one-dimensional case, the core script is
executed using:

``` r
library(doParallel)
nodes_CL = detectCores()   # Detect number of cores to use
cl = makeCluster(nodes_CL) # Specify number of threads here
registerDoParallel(cl)

source("1D_simulations/simulations_1D.R", local = TRUE)

stopCluster(cl)
```

The simulations can be executed in series by commenting the lines with
`doParallel` commands and sourcing the script as it is:

``` r
# library(doParallel)
# nodes_CL = detectCores()   
# cl = makeCluster(nodes_CL) 
# registerDoParallel(cl)

source("1D_simulations/simulations_1D.R", local = TRUE)

# stopCluster(cl)
```
