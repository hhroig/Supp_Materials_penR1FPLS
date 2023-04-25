library(tidyverse)
library(gridExtra)
library(viridis)
library(ggpubr)
library(scales)


compare_methods_fun <- function(input_folder){

  out_folder <- paste0(input_folder, "results_plots/")

  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }

  # ALL RES files ---------------------------------------------------

  bind_all_res <- data.frame()

  for (f_in in list.files(path = input_folder, pattern = "all_final_mods_res")) {
    bind_all_res <- rbind(
      bind_all_res,
      readRDS(paste0(input_folder, f_in) )
    )
  }

  # bind_all_res <- readRDS(paste0(input_folder,"all_final_mods_res.Rds"))

  bind_all_res <- bind_all_res %>%
    as_tibble() %>%
    mutate(method = factor(method, levels = c("B-FPLS",
                                              "PB-FPLS",
                                              "MV-PLS",
                                              "R1-FPLS")),
           nComp = as.factor(nComp),
           Rsq = paste0("R^2 = ", Rsq),
           beta_num = paste0("Beta_", beta_num))

  imse_limits <- range(bind_all_res$IMSE)
  mse_limits <- range(bind_all_res$mse_yy)
  mse_cl_limits <- range(bind_all_res$mse_cl)
  cor_limits <- range(bind_all_res$cor_yy, bind_all_res$cor_cl)


  # CVEs files -----------------------------------------------------------------

  all_cve <- data.frame()

  for (f_in in list.files(path = input_folder, pattern = "all_cv_res")) {
    all_cve <- rbind(
      all_cve,
      readRDS(paste0(input_folder, f_in) )
    )
  }

  # all_cve <- readRDS(paste0(input_folder, "all_cv_res.Rds"))

  all_cve <- all_cve %>%
    as_tibble() %>%
    mutate(method = factor(method, levels = c("B-FPLS",
                                              "PB-FPLS",
                                              "MV-PLS",
                                              "R1-FPLS")),
           nComp = as.factor(ncomp),
           Rsq = paste0("R^2 = ", Rsq),
           beta_num = paste0("Beta_", beta_num))

  # IMSE + MSE --------------------------------------------------------------

  color_codes <- c("B-FPLS" = hue_pal()(5)[2],
                   "PB-FPLS" = hue_pal()(5)[3],
                   "MV-PLS" = hue_pal()(5)[4],
                   "R1-FPLS" = hue_pal()(5)[1])

  for (rsq in unique(bind_all_res$Rsq)) {
    for (beta.num in unique(bind_all_res$beta_num)) {

      p_cve <- ggplot(all_cve %>% filter(Rsq == rsq, beta_num == beta.num),
                      aes(x = nComp, y = CVE, fill = method)) +
        geom_boxplot(position=position_dodge(0.8))  +
        ylab("CVE") +
        xlab("# of components") +
        scale_fill_manual(values = color_codes)+
        theme_bw()  +
        theme(legend.position="bottom", text = element_text(size = 20)) +
        labs(fill = "")

      p_imse <- ggplot(bind_all_res %>% filter(Rsq == rsq, beta_num == beta.num,
                                               method != "MV-PLS"),
                       aes(x = nComp, y = IMSE, fill = method)) +
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
                               paste0("2d_cve_imse_", rsq, beta.num,".png")  ),
             width = 12, height = 6 )
      ggsave(p_both,
             filename = paste0(out_folder,
                               paste0("2d_cve_imse_", rsq, beta.num,".pdf")  ),
             width = 12, height = 6 )



      # Log-scale
      p_cve_log <- ggplot(all_cve %>% filter(Rsq == rsq, beta_num == beta.num),
                          aes(x = nComp, y = log(CVE, base = 10), fill = method)) +
        geom_boxplot(position=position_dodge(0.8))  +
        ylab("log(CVE)") +
        xlab("# of components") +
        scale_fill_manual(values = color_codes)+
        theme_bw()  +
        theme(legend.position="bottom", text = element_text(size = 20)) +
        labs(fill = "")

      p_imse_log <- ggplot(bind_all_res %>% filter(Rsq == rsq, beta_num == beta.num,
                                                   method != "MV-PLS"),
                           aes(x = nComp, y = log(IMSE, base = 10), fill = method)) +
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
                               paste0("2d_cve_imse_log_", rsq, beta.num,
                                      ".pdf")  ),
             width = 12, height = 6 )
      ggsave(p_both_log,
             filename = paste0(out_folder,
                               paste0("2d_cve_imse_log_",  rsq, beta.num,
                                      ".png")  ),
             width = 12, height = 6 )



    } # loop "beta.num"
  }# loop "rsq"

  # Penalties ---------------------------------------------------------------

  df <- all_cve %>%
    filter(method %in% c("PB-FPLS",
                         "R1-FPLS"))

  for (beta.num in unique(df$beta_num)) {

    p_lamb <- ggplot(df %>% filter(beta_num == beta.num),
                     aes(x = nComp, y = log(penalty, base = 10))) +
      facet_grid(method~Rsq) +
      geom_boxplot() +
      ylab(expression(log(lambda ))) +
      xlab("# of components") +
      theme_bw() +
      theme(text = element_text(size = 20))

    # p_lamb

    ggsave(p_lamb,
           filename = paste0(out_folder, "2d_best_lambdas_",  beta.num, ".png"),
           width = 15, height = 10 )

    ggsave(p_lamb,
           filename = paste0(out_folder, "2d_best_lambdas_",  beta.num, ".pdf"),
           width = 15, height = 10 )
  }




  # Betas -------------------------------------------------------------------

  bind_all_betas <- data.frame()

  for (f_in in list.files(path = input_folder, pattern = "all_beta_hats_res")) {
    bind_all_betas <- rbind(
      bind_all_betas,
      readRDS(paste0(input_folder, f_in) )
    )
  }


  summ_all_betas <- bind_all_betas %>%
    as_tibble() %>%
    mutate(method = factor(method, levels = c("B-FPLS",
                                              "PB-FPLS",
                                              "MV-PLS",
                                              "R1-FPLS",
                                              "True Beta") )) %>%
    group_by(method, Rsq, beta_num, nComp, x, y) %>%
    summarise(mean_z = mean(z))


  plot_mean_betas <- function(summ_all_betas,
                              beta.num = 3,
                              n.Comp = 4,
                              rsq.in = 0.9,
                              bins = NA,
                              bin.width = 0.5,
                              line.size = 1) {

    betas_limits <-  summ_all_betas %>%
      filter(beta_num == beta.num,
             Rsq == rsq.in) %>%
      ungroup() %>%
      select(mean_z) %>%
      range()

    plot_data <- summ_all_betas %>%
      filter(beta_num == beta.num, nComp == n.Comp)

    betas_plot <- ggplot(data = plot_data,
                         mapping = aes(x = x, y = y, z = mean_z)) +
      geom_raster(aes(fill = mean_z)) +
      facet_grid( ~ method) +
      scale_fill_viridis() +
      theme_bw()+
      coord_fixed()+
      labs(fill = "") + xlab("p1") + ylab("p2") +
      theme(text = element_text(size = 20))
    # print(betas_plot)

    betas_plot_cont <- ggplot(data = plot_data,
                              mapping = aes(x = x, y = y, z = mean_z)) +
      geom_contour(aes(color=after_stat(level)), linewidth = line.size,
                   bins = bins, binwidth = bin.width) +
      facet_grid( ~ method) +
      scale_color_viridis() +
      theme_bw() +
      coord_fixed()+xlab("p1") + ylab("p2") +
      theme(text = element_text(size = 20))
    # print(betas_plot_cont)

    return(list(
      p_both_fill = betas_plot,
      p_both_cont = betas_plot_cont
    ))

  }

  out_folder_mean_betas <- paste0(out_folder, "mean_betas_no_mv/")

  if (!dir.exists(out_folder_mean_betas)) {
    dir.create(out_folder_mean_betas)
  }


  for (rsq.in in unique(summ_all_betas$Rsq) ) {

    for (n.Comp in unique(summ_all_betas$nComp)) {

      for (n.Beta in unique(summ_all_betas$beta_num)) {


        p_mean_fill <- plot_mean_betas(summ_all_betas %>% filter(method != "MV-PLS"),
                                       beta.num = n.Beta,
                                       n.Comp = n.Comp,
                                       rsq.in = rsq.in,
                                       bins = NULL,
                                       bin.width = 0.1,
                                       line.size = 0.8)[["p_both_fill"]]

        ggsave(p_mean_fill,
               filename = paste0(out_folder_mean_betas,
                                 "2d_mean_beta_", n.Beta,
                                 "_nComp_", n.Comp,
                                 "_Rsq_", rsq.in,
                                 ".png"),
               width = 12, height = 6 )
        ggsave(p_mean_fill,
               filename = paste0(out_folder_mean_betas,
                                 "2d_mean_beta_", n.Beta,
                                 "_nComp_", n.Comp,
                                 "_Rsq_", rsq.in,
                                 ".pdf"),
               width = 12, height = 6 )

      } # loop beta_num
    } # loop nComp
  } # loop rsq.in


  ## Just MV-PLS Betas ----

  plot_mean_mv_betas <- function(summ_all_betas,
                                 n.Comp = 4,
                                 rsq.in = 0.9,
                                 bins = NA,
                                 bin.width = 0.5,
                                 line.size = 1) {

    betas_limits <-  summ_all_betas %>%
      filter(Rsq == rsq.in) %>%
      ungroup() %>%
      select(mean_z) %>%
      range()

    plot_data <- summ_all_betas %>%
      filter( nComp == n.Comp)

    betas_plot <- ggplot(data = plot_data,
                         mapping = aes(x = x, y = y, z = mean_z)) +
      geom_raster(aes(fill = mean_z)) +
      facet_grid( method ~ beta_name) +
      scale_fill_viridis() +
      theme_bw()+
      coord_fixed()+
      labs(fill = "") + xlab("p1") + ylab("p2") +
      theme(text = element_text(size = 20))
    # print(betas_plot)


    return(list(
      p_both_fill = betas_plot
    ))

  }

  out_folder_mean_betas <- paste0(out_folder, "mean_betas_just_mv/")

  if (!dir.exists(out_folder_mean_betas)) {
    dir.create(out_folder_mean_betas)
  }

  summ_all_betas_mv <- summ_all_betas %>%
    filter(method == "MV-PLS") %>%
    mutate(beta_name = case_when(
      beta_num == 3 ~ "Double exponential",
      beta_num == 5 ~ "Monkey saddle",
      beta_num == 6 ~ "Rough",
      .default = as.character(beta_num)
    ) )

  for (rsq.in in unique(summ_all_betas$Rsq) ) {
    for (n.Comp in unique(summ_all_betas$nComp)) {


      p_mean_fill <- plot_mean_mv_betas(summ_all_betas_mv,
                                        n.Comp = n.Comp,
                                        rsq.in = rsq.in,
                                        bins = NULL,
                                        bin.width = 0.1,
                                        line.size = 0.8)[["p_both_fill"]]

      ggsave(p_mean_fill,
             filename = paste0(out_folder_mean_betas,
                               "2d_mean_mv_beta_",
                               "_nComp_", n.Comp,
                               "_Rsq_", rsq.in,
                               ".png"),
             width = 8, height = 6 )

      ggsave(p_mean_fill,
             filename = paste0(out_folder_mean_betas,
                               "2d_mean_mv_beta_",
                               "_nComp_", n.Comp,
                               "_Rsq_", rsq.in,
                               ".pdf"),
             width = 8, height = 6 )

    } # loop nComp
  } # loop rsq.in




}# end function
