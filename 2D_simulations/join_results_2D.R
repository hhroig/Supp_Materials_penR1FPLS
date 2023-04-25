all_files <- list.files(input_dir, pattern = "_rep_")

# Useful funs -----------------------------------------------------------------

mass_mat_fun <- fdaPDE:::CPP_get.FEM.Mass.Matrix

imse_beta_2d <- penR1FPLS:::imse_beta_2d

get_set_param <- function(string){

  rep_num <- str_sub(string,
                     start = str_locate(string, "_rep_")[ , 2]+1, # where _rep_ ends + 1
                     end = str_locate(string, ".Rds")[ , 1]-1  ) # where .Rds starts - 1

  beta_num <- str_sub(string,
                      start = str_locate(string, "_beta_")[ , 2]+1, # where _beta_ ends + 1
                      end = str_locate(string, "_rep_")[ , 1]-1  ) # where _rep_ starts - 1

  Rsq <- str_sub(string,
                 start = str_locate(string, "R2_")[ , 2]+1, # where _R2_ ends + 1
                 end = str_locate(string, "_beta_")[ , 1]-1  ) # where _beta_ starts - 1

  method <- str_sub(string,
                    start = 1, # from the start
                    end = str_locate(string, "_R2_")[ , 1]-1  ) # where _R2_ ends - 1

  return(
    list(
      rep_num = as.numeric(rep_num),
      beta_num = as.numeric(beta_num),
      Rsq = as.numeric(Rsq),
      method = method
    )
  )

}



# Loop in files -----------------------------------------------------------


load(paste0(input_dir, "settings.RData"))
R0tps <- readRDS(paste0(input_dir, "R0tps.Rds"))

total_reps <- max(get_set_param(all_files)$rep_num)
rep_starts <- min(get_set_param(all_files)$rep_num)

beta_num_lst <- unique(get_set_param(all_files)$beta_num)
Rsq_lst <- unique(get_set_param(all_files)$Rsq)


all_cv_results <- data.frame()
all_beta_hats <- data.frame()
all_res <- data.frame()



for (rep_num in rep_starts:total_reps) {
  for(Rsq in Rsq_lst){
    for (beta_num in beta_num_lst) {


      cat("-->  Joining R^2=" , Rsq, " [", Rsq_lst, "]",
          " | Beta ", beta_num, " [", beta_num_lst, "]",
          " | Rep ", rep_num, "/", total_reps, "\n")


      rep_ident <- paste0("_R2_", Rsq,
                          "_beta_", beta_num,
                          "_rep_", rep_num)

      ## Load data ------------------------------------------------------------

      L <- readRDS(paste0(input_dir, "data", rep_ident, ".Rds"))

      X_true <- L[["X"]]
      Y_true <- L[["Y"]]
      Y_clean <- L[["Y_clean"]] %>% as.numeric()
      beta_true <-  L[["coefficient_function"]] %>% as.numeric()

      # FE basis, mesh, and mass R0 (just to integrate in IMSE)
      FEM_basis <- L[["basisobj"]]
      mesh <- L[["mesh"]]
      nodes <- mesh$nodes

      R0 <- mass_mat_fun(FEMbasis = FEM_basis)

      ## CVEs ------------------------------------------------------------

      cv_R1FPLS <- readRDS(paste0(input_dir, "R1-FPLS", rep_ident, ".Rds"))
      cv_PBFPLS <- readRDS(paste0(input_dir, "PB-FPLS", rep_ident, ".Rds"))
      cv_BFPLS <- readRDS(paste0(input_dir, "B-FPLS", rep_ident, ".Rds"))
      cv_MV <- readRDS(paste0(input_dir, "MV-PLS", rep_ident, ".Rds"))


      res_cv_R1FPLS <- data.frame(data.frame(CVE = cv_R1FPLS[["CVEs_ncomp"]],
                                             ncomp = 1:length(cv_R1FPLS[["CVEs_ncomp"]])),
                                  data.frame(penalty = cv_R1FPLS[["best_penalties"]]),
                                  method = "R1-FPLS",
                                  Rsq = Rsq,
                                  beta_num = beta_num,
                                  rep_num = rep_num  )


      res_cv_PBFPLS <- data.frame(data.frame(CVE = cv_PBFPLS[["CVEs_ncomp"]],
                                             ncomp = 1:length(cv_PBFPLS[["CVEs_ncomp"]])),
                                  data.frame(penalty = cv_PBFPLS[["best_penalties"]]),
                                  method = "PB-FPLS",
                                  Rsq = Rsq,
                                  beta_num = beta_num,
                                  rep_num = rep_num  )

      res_cv_BFPLS <- data.frame(data.frame(CVE = cv_BFPLS[["CVEs_ncomp"]],
                                            ncomp = 1:length(cv_BFPLS[["CVEs_ncomp"]])),
                                 penalty = NA,
                                 method = "B-FPLS",
                                 Rsq = Rsq,
                                 beta_num = beta_num,
                                 rep_num = rep_num  )

      res_cv_MV <- data.frame(data.frame(CVE = cv_MV[["CVEs_ncomp"]],
                                         ncomp = 1:length(cv_MV[["CVEs_ncomp"]])),
                              penalty = NA,
                              method = "MV-PLS",
                              Rsq = Rsq,
                              beta_num = beta_num,
                              rep_num = rep_num  )

      # Save all CVES and best penalties into single matrix:
      all_cv_results <- rbind(
        all_cv_results,
        res_cv_R1FPLS ,
        res_cv_PBFPLS ,
        res_cv_BFPLS,
        res_cv_MV   )


      ## Loop in nComp -----------------------------------------------------------


      for (nComp in 1:max_nComp) {

        cat("          ncomp = ", nComp, "/", max_nComp, "\n")

        df_beta <- data.frame(x = nodes[, 1],
                              y = nodes[, 2],
                              z = beta_true,
                              method = "True Beta",
                              Rsq = Rsq,
                              beta_num = beta_num,
                              nComp = nComp,
                              rep_num = rep_num)
        all_beta_hats <- rbind(all_beta_hats, df_beta)

        ## FPLS -----------------------------------------------------------------

        final_FPLS <- cv_R1FPLS[["final_model"]]

        Y_hat <- fitted.values(final_FPLS)[, , nComp] %>% as.numeric()

        beta_hat <- final_FPLS[["coefficient_function"]][, , nComp] %>% as.numeric()

        df_beta <- data.frame(x = nodes[, 1],
                              y = nodes[, 2],
                              z = beta_hat,
                              method = "R1-FPLS",
                              Rsq = Rsq,
                              beta_num = beta_num,
                              nComp = nComp,
                              rep_num = rep_num)
        # Metrics:

        res_fpls <- data.frame(IMSE = imse_beta_2d(beta_hat,
                                                   beta_true, R0, 1),
                               cor_yy = cor(as.numeric(Y_true), Y_hat),
                               cor_cl = cor(Y_clean, Y_hat),
                               mse_yy = mean( (as.numeric(Y_true) - Y_hat)^2 ),
                               mse_cl = mean( (Y_clean - Y_hat)^2 ),
                               method = "R1-FPLS",
                               Rsq = Rsq,
                               beta_num = beta_num,
                               nComp = nComp,
                               rep_num = rep_num)
        # Bind FPLS results

        all_res <- rbind(all_res, res_fpls)
        all_beta_hats <- rbind(all_beta_hats, df_beta)

        ## pTPS ------------------------------------------------------------------

        ptps_sol <- cv_PBFPLS[["final_model"]]

        Y_hat <- fitted.values(ptps_sol)[, , nComp] %>% as.numeric()

        beta_hat <- ptps_sol[["coefficient_function"]][, , nComp] %>% as.numeric()

        df_beta <- data.frame(x = nodes[, 1],
                              y = nodes[, 2],
                              z = beta_hat,
                              method = "PB-FPLS",
                              Rsq = Rsq,
                              beta_num = beta_num,
                              nComp = nComp,
                              rep_num = rep_num)
        # Metrics:
        res_ptps <- data.frame(IMSE = imse_beta_2d(beta_hat,
                                                   beta_true, R0, 1),
                               cor_yy = cor(as.numeric(Y_true), Y_hat),
                               cor_cl = cor(Y_clean, Y_hat),
                               mse_yy = mean( (as.numeric(Y_true) - Y_hat)^2 ),
                               mse_cl = mean( (Y_clean - Y_hat)^2 ),
                               method = "PB-FPLS",
                               Rsq = Rsq,
                               beta_num = beta_num,
                               nComp = nComp,
                               rep_num = rep_num)
        # Bind FPLS results
        all_res <- rbind(all_res, res_ptps)
        all_beta_hats <- rbind(all_beta_hats, df_beta)

        ## TPS -------------------------------------------------------------------

        tps_sol <- cv_PBFPLS[["final_model"]]

        Y_hat <- fitted.values(tps_sol)[, , nComp] %>% as.numeric()

        beta_hat <- tps_sol[["coefficient_function"]][, , nComp] %>% as.numeric()

        df_beta <- data.frame(x = nodes[, 1],
                              y = nodes[, 2],
                              z = beta_hat,
                              method = "B-FPLS",
                              Rsq = Rsq,
                              beta_num = beta_num,
                              nComp = nComp,
                              rep_num = rep_num)
        # Metrics:
        res_tps <- data.frame(IMSE = imse_beta_2d(beta_hat,
                                                  beta_true, R0, 1),
                              cor_yy = cor(as.numeric(Y_true), Y_hat),
                              cor_cl = cor(Y_clean, Y_hat),
                              mse_yy = mean( (as.numeric(Y_true) - Y_hat)^2 ),
                              mse_cl = mean( (Y_clean - Y_hat)^2 ),
                              method = "B-FPLS",
                              Rsq = Rsq,
                              beta_num = beta_num,
                              nComp = nComp,
                              rep_num = rep_num)
        # Bind FPLS results
        all_res <- rbind(all_res, res_tps)
        all_beta_hats <- rbind(all_beta_hats, df_beta)

        # MV-PLS -----------------------------------------------------------------

        mv_sol <- pls::plsr(Y_true ~ X_true,
                            ncomp =  nComp,
                            model = "oscorespls",
                            center = TRUE )

        Y_hat <-  predict(mv_sol, newdata = X_true, ncomp = nComp) %>% as.numeric()

        beta_hat <- coef(mv_sol, ncomp = nComp) %>% as.numeric()

        df_beta <- data.frame(x = nodes[, 1],
                              y = nodes[, 2],
                              z = beta_hat,
                              method = "MV-PLS",
                              Rsq = Rsq,
                              beta_num = beta_num,
                              nComp = nComp,
                              rep_num = rep_num)
        # Metrics:
        res_mvpls <- data.frame(IMSE = imse_beta_2d(beta_hat,
                                                    beta_true, R0, 1),
                                cor_yy = cor(as.numeric(Y_true), Y_hat),
                                cor_cl = cor(Y_clean, Y_hat),
                                mse_yy = mean( (as.numeric(Y_true) - Y_hat)^2 ),
                                mse_cl = mean( (Y_clean - Y_hat)^2 ),
                                method = "MV-PLS",
                                Rsq = Rsq,
                                beta_num = beta_num,
                                nComp = nComp,
                                rep_num = rep_num)
        # Bind FPLS results
        all_res <- rbind(all_res, res_mvpls)
        all_beta_hats <- rbind(all_beta_hats, df_beta)

      }  # end nComp loop




    } # end beta_num loop
  } # end Rsq loop
} # end rep_num loop


saveRDS(all_cv_results, file = paste0(input_dir, "all_cv_res.Rds"))
saveRDS(all_res, file = paste0(input_dir, "all_final_mods_res.Rds"))
saveRDS(all_beta_hats, file = paste0(input_dir, "all_beta_hats_res.Rds"))
