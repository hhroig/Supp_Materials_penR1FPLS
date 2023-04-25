

# Tools -------------------------------------------------------------------

imse_beta_1d <- function(argvals, beta_true, beta_hat) {

  imse <- MESS::auc(argvals, (beta_true - beta_hat)^2,
                    type = "linear")/abs(diff(range(argvals)))

  return(imse)

}

# Fixed X data ------------------------------------------------------------

# # Octane number:
# Y <- as.matrix(gasoline$octane)

# Gasoline NIR spectra:
X <- as.matrix( as.data.frame(gasoline$NIR) )
Xc <- scale(X, center = TRUE, scale = FALSE)

# Wavenumber:
argvals <- seq(900, 1700, by = 2)

# Beta function -----------------------------------------------------------

beta_p <- function(argvals){
  t <- (argvals - 900)/800
  y <- 2*sin(0.5*pi*t) + 4*sin(1.5*pi*t) + 5*sin(2.5*pi*t)
  return(y)
}

beta_true <- beta_p(argvals)


# Generate Y  (clean) using numerical integration:
Y_clean <- matrix(NA, nrow = nrow(Xc), ncol = 1)

for (j in 1:nrow(Xc)) {
  Y_clean[j] <-  MESS::auc(argvals, Xc[j, ]*beta_true,
                           type = 'linear')
}




# Basis functions settings -----

# Ruppert's law:  nbasis = nbreaks + norder - 2  and norder = degree + 1
n_breaks <- min(round(length(argvals)/4), 40)
n_basis <- n_breaks + (3+1) - 2

# B-spline basis:
bs_basis <- create.bspline.basis(rangeval = range(argvals),
                                 nbasis = n_basis)


save(lambdas_in, total_reps, max_nComp, my_setting,
     file = paste0(out_folder, "sim_setting.RData"))


# Repetitions -----

for(row_setting in 1:nrow(my_setting))       {
  for(rep_num in rep_starts:total_reps)

  {

    Rsq <- my_setting[row_setting, "Rsq"]

    # Initialize savings:
    all_CVEs <- data.frame()
    all_beta_hats <- data.frame()
    all_final_res <- data.frame()
    best_lambdas <- data.frame()

    cat("Doing setting ", row_setting, "/", nrow(my_setting), "\n")
    cat("     rep ", rep_num, "/", total_reps, "\n")


    # Noisy Y -----------------------------------------------------------

    var_e <- (1/Rsq - 1)*var(Y_clean)   # variance of errors

    Y <- Y_clean +
      as.matrix(
        rnorm(length(Y_clean), mean = 0, sd = sqrt(var_e))
      )


    folds <- caret::createFolds(Y, k = num_folds)

    # CV R1-FPLS ------------------------------------------------------------------

    cv_r1fpls <- cv_unique_par(X = X, Y = Y, argvals = argvals, penalty_vec = lambdas_in,
                               ncomp = max_nComp, folds = folds, basisobj = bs_basis,
                               verbose = TRUE, stripped = FALSE, method = "r1fpls_bs")

    cat("      cv_r1fpls done\n")

    best_lambdas <- rbind(
      best_lambdas,
      data.frame(lambda = cv_r1fpls$best_penalties,
                 method = "R1-FPLS",
                 nComp = 1:max_nComp,
                 rep_num = rep_num,
                 Rsq = Rsq)
    )


    all_CVEs <- rbind(
      all_CVEs,
      data.frame(
        CVE = as.numeric(  cv_r1fpls$CVEs_ncomp  ),
        method = "R1-FPLS",
        nComp = 1:max_nComp,
        rep_num = rep_num,
        Rsq = Rsq
      )
    )


    # CV PB-FPLS ------------------------------------------------------------------

    cv_pbfpls <- cv_unique_par(X = X, Y = Y, argvals = argvals, penalty_vec = lambdas_in,
                               ncomp = max_nComp, folds = folds, basisobj = bs_basis,
                               verbose = TRUE, stripped = FALSE, method = "fpls_bs")
    cat("      cv_pbfpls done\n")

    best_lambdas <- rbind(
      best_lambdas,
      data.frame(lambda = as.numeric( cv_pbfpls$best_penalties ),
                 method = "PB-FPLS",
                 nComp = 1:max_nComp,
                 rep_num = rep_num,
                 Rsq = Rsq)
    )

    all_CVEs <- rbind(
      all_CVEs,
      data.frame(
        CVE = as.numeric( cv_pbfpls$CVEs_ncomp ),
        method = "PB-FPLS",
        nComp = 1:max_nComp,
        rep_num = rep_num,
        Rsq = Rsq
      )
    )

    # CV B-FPLS (no penalties) ---------------------------------------------------

    cv_bfpls <- cv_unique_par(X = X, Y = Y, argvals = argvals, penalty_vec = 0,
                              ncomp = max_nComp, folds = folds, basisobj = bs_basis,
                              verbose = TRUE, stripped = FALSE, method = "fpls_bs")
    cat("      cv_bfpls done\n")

    best_lambdas <- rbind(
      best_lambdas,
      data.frame(lambda = as.numeric(cv_bfpls$best_penalties  ),
                 method = "B-FPLS",
                 nComp = 1:max_nComp,
                 rep_num = rep_num,
                 Rsq = Rsq)
    )

    all_CVEs <- rbind(
      all_CVEs,
      data.frame(
        CVE = as.numeric(cv_bfpls$CVEs_ncomp  ),
        method = "B-FPLS",
        nComp = 1:max_nComp,
        rep_num = rep_num,
        Rsq = Rsq
      )
    )


    # Final models for each ncomp ---------------------------------------------


    for (nComp in 1:max_nComp) {
      cat("       final models cmpt:", nComp, "/", max_nComp, "\n")

      ## R1-FPLS model ----

      m_final <- cv_r1fpls[["final_model"]]

      all_beta_hats <- rbind(all_beta_hats,
                             data.frame(argvals = argvals,
                                        z = m_final$coefficient_function[, , nComp],
                                        z_true = beta_true,
                                        method = "R1-FPLS",
                                        Rsq = Rsq,
                                        nComp = nComp,
                                        rep_num = rep_num) )

      all_final_res <- rbind(all_final_res,
                             tibble(
                               imse = imse_beta_1d(argvals = argvals,
                                                   beta_true = beta_true,
                                                   beta_hat = m_final$coefficient_function[, , nComp]),
                               correl = cor(as.numeric(Y), as.numeric(m_final$fitted.values[, , nComp])),
                               correl_clean = cor(as.numeric(Y_clean), as.numeric(m_final$fitted.values[, , nComp])),
                               correl2 = correl^2,
                               correl2_clean = correl_clean^2,
                               MSE = mean( as.numeric(Y - m_final$fitted.values[, , nComp])^2 ),
                               MSE_clean = mean( as.numeric(Y_clean - m_final$fitted.values[, , nComp])^2 ),
                               Rsq = Rsq,
                               nComp = nComp,
                               rep_num = rep_num,
                               method = "R1-FPLS")
      )

      rm(m_final) # I'll reuse the same model name

      ## PB-FPLS model ----


      m_final <- cv_pbfpls[["final_model"]]

      all_beta_hats <- rbind(all_beta_hats,
                             data.frame(argvals = argvals,
                                        z = as.numeric(m_final$coefficient_function[, , nComp]  ),
                                        z_true = beta_true,
                                        method = "PB-FPLS",
                                        Rsq = Rsq,
                                        nComp = nComp,
                                        rep_num = rep_num) )

      all_final_res <- rbind(all_final_res,
                             tibble(
                               imse = imse_beta_1d(argvals = argvals,
                                                   beta_true = beta_true,
                                                   beta_hat = m_final$coefficient_function[, , nComp]),
                               correl = cor(as.numeric(Y), as.numeric(m_final$fitted.values[, , nComp])),
                               correl_clean = cor(as.numeric(Y_clean), as.numeric(m_final$fitted.values[, , nComp])),
                               correl2 = correl^2,
                               correl2_clean = correl_clean^2,
                               MSE = mean( as.numeric(Y - m_final$fitted.values[, , nComp])^2 ),
                               MSE_clean = mean( as.numeric(Y_clean - m_final$fitted.values[, , nComp])^2 ),
                               Rsq = Rsq,
                               nComp = nComp,
                               rep_num = rep_num,
                               method = "PB-FPLS")
      )

      rm(m_final) # I'll reuse the same model name



      ## B-FPLS model ----

      m_final <- cv_bfpls[["final_model"]]

      all_beta_hats <- rbind(all_beta_hats,
                             data.frame(argvals = argvals,
                                        z = as.numeric(m_final$coefficient_function[, , nComp]  ),
                                        z_true = beta_true,
                                        method = "B-FPLS",
                                        Rsq = Rsq,
                                        nComp = nComp,
                                        rep_num = rep_num) )

      all_final_res <- rbind(all_final_res,
                             tibble(
                               imse = imse_beta_1d(argvals = argvals,
                                                   beta_true = beta_true,
                                                   beta_hat = m_final$coefficient_function[, , nComp]),
                               correl = cor(as.numeric(Y), as.numeric(m_final$fitted.values[, , nComp])),
                               correl_clean = cor(as.numeric(Y_clean), as.numeric(m_final$fitted.values[, , nComp])),
                               correl2 = correl^2,
                               correl2_clean = correl_clean^2,
                               MSE = mean( as.numeric(Y - m_final$fitted.values[, , nComp])^2 ),
                               MSE_clean = mean( as.numeric(Y_clean - m_final$fitted.values[, , nComp])^2 ),
                               Rsq = Rsq,
                               nComp = nComp,
                               rep_num = rep_num,
                               method = "B-FPLS")
      )

      rm(m_final) # I'll reuse the same model name


    }  # end nComp loop


    # Savings:

    saveRDS(all_CVEs, file =
              paste0(out_folder,
                     "cves_set_",
                     row_setting,
                     "_rep_",
                     rep_num,
                     ".Rds"))

    saveRDS(all_final_res, file =
              paste0(out_folder,
                     "final_models_set_",
                     row_setting,
                     "_rep_",
                     rep_num,
                     ".Rds"))

    saveRDS(all_beta_hats, file =
              paste0(out_folder,
                     "betas_set_",
                     row_setting,
                     "_rep_",
                     rep_num,
                     ".Rds"))

    saveRDS(best_lambdas, file =
              paste0(out_folder,
                     "best_lambdas_set_",
                     row_setting,
                     "_rep_",
                     rep_num,
                     ".Rds"))

  } # end nested inner
} # end nested outer


