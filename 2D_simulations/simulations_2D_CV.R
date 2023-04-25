

# Useful funs -----------------------------------------------------------------
wlog = function(text,...){

  cat(paste0(date(),"	", text,...,"
"), file = paste0(out_folder, "log.txt"), append = T)

  cat( paste0(date(),"	", text, ...) )
}


# Initial savings: --------------------------------------------------------

save(lambdas_in,
     num_of_samples,
     n_basis_tps,
     total_reps,
     max_nComp,
     Rsq_lst,
     beta_num_lst,
     file = paste0(out_folder, "settings.RData"))


# Fix R0 TPS --------------------------------------------------------------

L <- generate_2d_data(x, y,
                      num_samples = 2,
                      beta_num = 3,
                      Rsq = 1)
X_true <- L[["X"]]
mesh <- L[["mesh"]]
nodes <- mesh$nodes

gam_fit <- mgcv::gam(X_true[1, ] ~ s(nodes[ , 1],
                                     nodes[ , 2],
                                     bs = "tp",
                                     k = n_basis_tps))
# Evaluate basis functions:
# (rows corresponding to argument values and columns to basis functions)
# X is approximated by A %*% t(Psi):
Psi <- stats::model.matrix(gam_fit)
tPsi <- Matrix::t(Psi)

# Matrix of inner products (mass):
R0tps <- matrix(NA, nrow = ncol(Psi), ncol = ncol(Psi))

# numerical approx. of the inner products:
for (i in 1:nrow(R0tps)) {
  for (j in i:ncol(R0tps)) {
    df <- as.data.frame(nodes)
    df$z = as.numeric(Psi[, i]*Psi[, j])
    R0tps[i,j] <-  penR1FPLS:::getVolume(df)
  } #j
} #i

R0tps[lower.tri(R0tps)] <- R0tps[upper.tri(R0tps)]

saveRDS(R0tps, file = paste0(out_folder, "R0tps.Rds"))

rm(gam_fit, mesh, nodes, L, X_true)


# Repetitions -------------------------------------------------------------

for (rep_num in rep_starts:total_reps) {


  for(Rsq in Rsq_lst){
    for (beta_num in beta_num_lst) {

      rep_ident <- paste0("_R2_", Rsq,
                          "_beta_", beta_num,
                          "_rep_", rep_num)

      wlog("R^2=" , Rsq, " [", Rsq_lst, "]",
           " | Beta ", beta_num, " [", beta_num_lst, "]",
           " | Rep ", rep_num, "/", total_reps, "\n")

      ## Generate data ---------------------------------------------------------

      L <- generate_2d_data(x, y,
                            num_samples = num_of_samples,
                            beta_num = beta_num,
                            Rsq = Rsq)

      X_true <- L[["X"]]
      Y_true <- L[["Y"]]
      Y_clean <- L[["Y_clean"]] %>% as.numeric()
      beta_true <-  L[["coefficient_function"]] %>% as.numeric()

      # FE basis, mesh, and mass R0 (just to integrate in IMSE)
      FEM_basis <- L[["basisobj"]]
      mesh <- L[["mesh"]]
      nodes <- mesh$nodes


      # Same CVE folds for all:

      folds <- caret::createFolds(Y_true, k = k_folds)

      saveRDS(L, file = paste0(out_folder, "data", rep_ident, ".Rds"))



      ## R1FPLS ----------------------------------------------------------------

      wlog("      R1-FPLS\n")

      cv_r1fpls <- cv_unique_par(X = X_true,
                                 Y = Y_true,
                                 penalty_vec = lambdas_in,
                                 ncomp = max_nComp,
                                 folds = folds,
                                 basisobj = FEM_basis,
                                 method = "r1fpls_fem",
                                 verbose = FALSE,
                                 stripped = FALSE )

      saveRDS(cv_r1fpls, file = paste0(out_folder, "R1-FPLS", rep_ident, ".Rds"))

      # Penalized TPS ---------------------------------------------------------

      wlog("      PB-FPLS\n")

      cv_pTPS <- cv_unique_par(X = X_true,
                               Y = Y_true,
                               center = TRUE,
                               nodes = nodes,
                               nbasis = n_basis_tps,
                               penalty_vec = lambdas_in,
                               basisobj = NA,
                               ncomp = max_nComp,
                               folds = folds,
                               method = "fpls_tps",
                               verbose = FALSE,
                               stripped = FALSE,
                               R0 = R0tps )
      saveRDS(cv_pTPS, file = paste0(out_folder, "PB-FPLS", rep_ident, ".Rds"))


      ## TPS no penalties ------------------------------------------------------

      wlog("      B-FPLS\n")

      cv_0_TPS <- cv_unique_par(X = X_true,
                                Y = Y_true,
                                center = TRUE,
                                nodes = nodes,
                                nbasis = n_basis_tps,
                                basisobj = NA,
                                penalty_vec = 0,
                                ncomp = max_nComp,
                                folds = folds,
                                method = "fpls_tps",
                                verbose = FALSE,
                                stripped = FALSE,
                                R0 = R0tps  )

      saveRDS(cv_0_TPS, file = paste0(out_folder, "B-FPLS", rep_ident, ".Rds"))

      # MV-PLS -----------------------------------------------------------------

      wlog("      MV-PLS\n")

      cv_mv <- cv_mvpls(X = X_true,
                        Y = Y_true,
                        ncomp = max_nComp,
                        center = TRUE,
                        folds = folds,
                        verbose = FALSE )

      saveRDS(cv_mv, file = paste0(out_folder, "MV-PLS", rep_ident, ".Rds"))

    } # end beta_num loop

  } # end Rsq loop

} # end rep_num loop


