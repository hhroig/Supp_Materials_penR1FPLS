library(tidyverse)
library(gridExtra)
library(ggpubr)
library(scales)



compare_methods_fun <- function(input_folder){

  out_folder <- paste0(input_folder, "results_plots_v2/")

  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }

  ## Data ---------------------------------------------------

  # Final Models Files

  all_final_res <- data.frame()

  final_res_files <- list.files(path = input_folder, pattern = "final_models")

  for (ind_file in final_res_files) {

    all_final_res <- rbind(
      all_final_res,
      readRDS(paste0(input_folder, ind_file))
    )

  }


  # CVEs files:

  all_cves <- data.frame()

  all_cves_files <- list.files(path = input_folder, pattern = "cves_")

  for (ind_file in all_cves_files) {

    all_cves <- rbind(
      all_cves,
      readRDS(paste0(input_folder, ind_file))
    )

  }


  # Best Lambdas files:

  all_best_lambdas <- data.frame()

  all_best_lambdas_files <- list.files(path = input_folder, pattern = "best_lambdas")

  for (ind_file in all_best_lambdas_files) {

    all_best_lambdas <- rbind(
      all_best_lambdas,
      readRDS(paste0(input_folder, ind_file))
    )

  }

  # Betas files:

  all_betas <- data.frame()

  all_betas_files <- list.files(path = input_folder, pattern = "betas_")

  for (ind_file in all_betas_files) {

    all_betas <- rbind(
      all_betas,
      readRDS(paste0(input_folder, ind_file))
    )

  }


  all_final_res <- all_final_res %>%
    mutate(nComp = as.factor(nComp))

  all_cves <- all_cves %>%
    mutate(nComp = as.factor(nComp))

  all_best_lambdas <- all_best_lambdas %>%
    mutate(nComp = as.factor(nComp))


  # Limits:

  mse_limits <- range(all_final_res$MSE, all_final_res$MSE_clean)
  cve_limits <- range(all_cves$CVE)
  corr_limits <- range(all_final_res$correl, all_final_res$correl_clean )
  imse_limits <- range(all_final_res$imse )


  # IMSE + MSE --------------------------------------------------------------


  color_codes <- c("B-FPLS" = hue_pal()(5)[2],
                   "PB-FPLS" = hue_pal()(5)[3],
                   # "R1-FPLS-FE" = hue_pal()(5)[5],
                   "R1-FPLS" = hue_pal()(5)[1])


  for (rsq in unique(all_final_res$Rsq)) {

    p_cve <- ggplot(all_cves %>% filter(Rsq == rsq),
                    aes(x = nComp, y = CVE, fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("CVE") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")

    p_imse <- ggplot(all_final_res %>% filter(Rsq == rsq),
                     aes(x = nComp, y = imse, fill = method)) +
      geom_boxplot(position=position_dodge(0.8)) +
      ylab("IMSE") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="none", text = element_text(size = 20))

    # get legend
    leg <- get_legend(p_cve)
    # remove from last plot
    p_cve <- p_cve + theme(legend.position="none")

    p_both <- grid.arrange(p_cve, p_imse, nrow = 1, bottom = leg)

    ggsave(p_both,
           filename = paste0(out_folder,
                             paste0("cve_imse_rsq_", rsq,".png")  ),
           width = 12, height = 6 )
    ggsave(p_both,
           filename = paste0(out_folder,
                             paste0("cve_imse_rsq_", rsq,".pdf")  ),
           width = 12, height = 6 )


    # Log-scale
    p_cve_log <- ggplot(all_cves %>% filter(Rsq == rsq),
                        aes(x = nComp, y = log(CVE, base = 10), fill = method)) +
      geom_boxplot(position=position_dodge(0.8))  +
      ylab("log(CVE)") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes)+
      theme_bw()  +
      theme(legend.position="bottom", text = element_text(size = 20)) +
      labs(fill = "")

    p_imse_log <- ggplot(all_final_res %>% filter(Rsq == rsq),
                         aes(x = nComp, y = log(imse, base = 10), fill = method)) +
      geom_boxplot( position=position_dodge(0.8) ) +
      ylab("log(IMSE)") +
      xlab("# of components") +
      scale_fill_manual(values = color_codes) +
      theme_bw() +
      theme(legend.position="none", text = element_text(size = 20))

    # get legend
    leg <- get_legend(p_cve_log)
    # remove from last plot
    p_cve_log <- p_cve_log + theme(legend.position="none")


    p_both_log <- grid.arrange(p_cve_log, p_imse_log, nrow = 1, bottom = leg)

    ggsave(p_both_log,
           filename = paste0(out_folder,
                             paste0("cve_imse_log_rsq_", rsq,
                                    ".pdf")  ),
           width = 12, height = 6 )
    ggsave(p_both_log,
           filename = paste0(out_folder,
                             paste0("cve_imse_log_rsq_", rsq,
                                    ".png")  ),
           width = 12, height = 6 )

  }# loop "rsq"



  # Betas -------------------------------------------------------------------

  summ_all_betas <- all_betas %>%
    as_tibble() %>%
    group_by(method, Rsq, nComp, argvals) %>%
    summarise(mean_z = mean(z),
              sd_z = sd(z),
              z_true = mean(z_true))


  plot_mean_betas <- function(summ_all_betas,
                              n.Comp = 2,
                              R.sq = 0.9) {

    betas_limits <-  summ_all_betas %>%
      ungroup() %>%
      dplyr::select(mean_z, z_true, sd_z) %>%
      range()

    plot_data <- summ_all_betas %>%
      filter(nComp == n.Comp, Rsq == R.sq)

    betas_plot <- ggplot(data = plot_data,
                         mapping = aes(x = argvals )) +
      geom_line(aes(y = z_true), color = "red", alpha = 0.5, linewidth = 0.8) +
      geom_line(aes(y = mean_z), color = "gray", alpha = 0.8, linewidth = 0.8) +
      geom_line(aes(y = mean_z + 2*sd_z), color = "gray",
                alpha = 0.8, linewidth = 0.8, linetype = "dashed") +
      geom_line(aes(y = mean_z - 2*sd_z), color = "gray",
                alpha = 0.8, linewidth = 0.8, linetype = "dashed") +
      facet_grid( ~ method) +
      labs(x = "p",
           y = expression(beta(p))) +
      theme_bw() +
      theme(text = element_text(size = 20))
    # print(betas_plot)

    betas_plot2 <- ggplot(data = plot_data,
                          mapping = aes(x = argvals )) +
      geom_line(aes(y = z_true), color = "red", alpha = 0.5, linewidth = 0.8) +
      geom_line(aes(y = mean_z), color = "gray", alpha = 0.8, linewidth = 0.8) +
      facet_grid( ~ method) +
      xlab("p") +
      ylab(expression(beta(p))) +
      theme_bw()
    # print(betas_plot2)

    return(list(betas_plot = betas_plot,
                betas_plot2 = betas_plot2))

  }



  plot_all_betas <- function(all_betas,
                             n.Comp = 2) {

    betas_limits <-  all_betas %>%
      ungroup() %>%
      dplyr::select(z, z_true) %>%
      range()

    plot_data <- all_betas %>%
      filter(nComp == n.Comp)

    betas_plot <- ggplot(data = plot_data,
                         mapping = aes(x = argvals, group = rep_num )) +
      geom_line(aes(y = z_true), color = "red", alpha = 0.8, linewidth = 0.8) +
      geom_line(aes(y = z), color = "gray", alpha = 0.5, linewidth = 0.8) +
      facet_grid( method ~ Rsq) +
      xlab("p") +
      ylab(expression(beta(p))) +
      theme_bw()
    # print(betas_plot)


    return(list(betas_plot = betas_plot))

  }

  ### Plot and save Betas ------------------------------------------------------

  out_folder_mean_betas <- paste0(out_folder, "mean_betas/")

  if (!dir.exists(out_folder_mean_betas)) {
    dir.create(out_folder_mean_betas)
  }



  for (n.Comp in unique(summ_all_betas$nComp)) {
    for (R.sq in unique(summ_all_betas$Rsq)) {

      p_mean <- plot_mean_betas(summ_all_betas,
                                n.Comp = n.Comp)[["betas_plot"]]

      # print(p_mean)
      ggsave(p_mean,
             filename = paste0(out_folder_mean_betas,
                               "mean_beta_nComp_", n.Comp,
                               "_Rsq_", R.sq,
                               ".png"),
             width = 12, height = 5 )
      ggsave(p_mean,
             filename = paste0(out_folder_mean_betas,
                               "mean_beta_nComp_", n.Comp,
                               "_Rsq_", R.sq,
                               ".pdf"),
             width = 12, height = 5 )

    }

  } # loop nComp plot Betas





  out_folder_all_betas <- paste0(out_folder, "all_betas/")

  if (!dir.exists(out_folder_all_betas)) {
    dir.create(out_folder_all_betas)
  }

  for (n.Comp in unique(all_betas$nComp)) {


    p_all <- plot_all_betas(all_betas,
                            n.Comp = n.Comp)[["betas_plot"]]

    # print(p_mean)
    ggsave(p_all,
           filename = paste0(out_folder_all_betas,
                             "all_betas_nComp_", n.Comp,
                             ".png"),
           width = 10, height = 10 )

  } # loop nComp plot Betas


  # Best lambdas ------------------------------------------------------------

  lambdas_limits <- range(all_best_lambdas$lambda)

  df <- all_best_lambdas %>%
    filter(method %in% c("PB-FPLS", "R1-FPLS")) %>%
    mutate(rep_num = as.factor(rep_num),
           nComp = as.factor(nComp))

  p_lamb <- ggplot(df, aes(x = nComp, y = log(lambda, base = 10))) +
    facet_grid(method~Rsq) +
    geom_boxplot() +
    ylab(expression(log(lambda ))) +
    xlab("# of components") +
    theme_bw() +
    theme(text = element_text(size = 20))

  ggsave(p_lamb,
         filename = paste0(out_folder, "best_lambdas_1D_sim.png"),
         width = 15, height = 10 )

  ggsave(p_lamb,
         filename = paste0(out_folder, "best_lambdas_1D_sim.pdf"),
         width = 15, height = 10 )

}# end function
