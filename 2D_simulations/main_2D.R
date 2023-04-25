# Please install the R package "penR1FPLS" from GitHub:
# devtools::install_github("hhroig/penR1FPLS", dependencies = TRUE)

library(dplyr)
library(tidyr)
library(fdaPDE)
library(penR1FPLS)
library(stringr)

# Settings ----------------------------------------------------------------


# length of the penalties grid:
num_lambdas <- 100

lambdas_in  <-  seq(-10, 12, length.out = num_lambdas)
lambdas_in  <-  10^(lambdas_in)

# number of TPS basis to consider:
n_basis_tps <- 10

# number of repetitions (total_reps - rep_starts)
rep_starts <- 1
total_reps  <-  100

# number K of folds to do cross-validation:
k_folds <- 5

# number of PLS components to compute:
max_nComp <- 10

# number of samples to generate:
num_of_samples  <-  60

# 2D domain grid:
num_grid_y_axes <- num_grid_x_axes <- 20
x <- seq(0, 1, length.out = num_grid_x_axes)
y <- seq(0, 1, length.out = num_grid_y_axes)

# Values of R^2 to consider (more than one is possible):
Rsq_lst <- 0.9

# True Betas (3 is "Double Exponential", 5 is "Monkey Saddle")
beta_num_lst  <-  c(3, 5)

# output folder:
out_folder <- "2D_simulations/results_2D/"
if (!dir.exists(out_folder)) {
  dir.create(out_folder)
}



# Call the actual simulations ----------------------------------------------

# Just the cross-validation:

library(doParallel)
nodes_CL = detectCores()   # Detect number of cores to use
cl = makeCluster(nodes_CL) # Specify number of threads here
registerDoParallel(cl)

source("2D_simulations/simulations_2D_CV.R", local = TRUE)

stopCluster(cl)

# Join results and compute summaries:
input_dir <- "2D_simulations/results_2D/"
source("2D_simulations/join_results_2D.R", local = TRUE)


# Plot comparisons --------------------------------------------------------


source("2D_simulations/compare_methods_2D.R", local = TRUE)

compare_methods_fun(input_folder = "2D_simulations/results_2D/")
