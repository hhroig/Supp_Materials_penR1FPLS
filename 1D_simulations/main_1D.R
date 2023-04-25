# Please install the R package "penR1FPLS" from GitHub:
# devtools::install_github("hhroig/penR1FPLS", dependencies = TRUE)

library(pls)
library(dplyr)
library(penR1FPLS)
library(fda)
library(fdaPDE)

# Settings ----------------------------------------------------------------

# length of the penalties grid:
num_lambdas <- 100

lambdas_in  <-  seq(-6, 12, length.out = num_lambdas)
lambdas_in  <-  10^(lambdas_in)

# number of repetitions (total_reps - rep_starts)
total_reps  <-  100
rep_starts <- 1

# number of PLS components to compute:
max_nComp <- 6

# number K of folds to do cross-validation:
num_folds <- 5

# different R^2 to fix (more than one are allowed):
my_setting <- expand.grid(
  Rsq = c(0.9)
)

# output folder:
out_folder <- "1D_simulations/results_1D/"
if (!dir.exists(out_folder)) {
  dir.create(out_folder)
}


# Call the actual simulations ----------------------------------------------

library(doParallel)
nodes_CL = detectCores()   # Detect number of cores to use
cl = makeCluster(nodes_CL) # Specify number of threads here
registerDoParallel(cl)

source("1D_simulations/simulations_1D.R", local = TRUE)

stopCluster(cl)


# Plot comparisons --------------------------------------------------------


source("1D_simulations/compare_methods_1D.R", local = TRUE)

compare_methods_fun(input_folder = "1D_simulations/results_1D/")
