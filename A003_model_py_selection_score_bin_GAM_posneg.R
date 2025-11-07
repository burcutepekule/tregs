rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(readr)  # For read_csv
library(stringr)
library(zoo)
library(mgcv)

source("/Users/burcutepekule/Dropbox/Treg_problem_v2/MISC/PLOT_FUNCTIONS.R")

## All saved already
# for(ss_start_threshold in c(50,100,150,200,250,300,350,400,450)){
#   source("~/Desktop/tregs/A002_analyse_py_selection_score_cts.R")
# }
sign_treshold = 0.001
df_params     = read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', show_col_types = FALSE)

for(inj_type in c('sterile','pathogenic')){
  # for(ss_start_threshold in c(50,100,150,200,250,300,350,400,450)){
  for(ss_start_threshold in c(250)){
    message("Processing ss_start_threshold ", ss_start_threshold, " inj type ",inj_type)
    source("~/Desktop/tregs/A002_analyse_py_selection_score_cts.R")
    df_raw     = readRDS(paste0('/Users/burcutepekule/Desktop/tregs/df_summary_selection_score_raw_',ss_start_threshold,'.rds'))
    
    df_model    = df_raw %>% dplyr::filter(comparison=='Treg_OFF_ON' & injury_type==inj_type)
    df_model    = df_model %>% dplyr::mutate(nonzero=ifelse(effect_size %in% c('Medium','Large') & abs(mean_diff)>tol, 1, 0))
    # df_model    = df_model %>% dplyr::mutate(nonzero=ifelse(abs(mean_diff)>tol, 1, 0))
    # df_model    = df_model %>% dplyr::mutate(nonzero=ifelse((mean_diff)>tol, 1, 0))
    
    df_model_e_neg = df_model %>% dplyr::filter(mean_diff<(-1)*tol)
    df_model_e_zer = df_model %>% dplyr::filter(mean_diff>=(-1)*tol & mean_diff<=(+1)*tol)
    df_model_e_pos = df_model %>% dplyr::filter(mean_diff>(+1)*tol)
    
    df_model    = merge(df_model, df_params, by='param_set_id')
    
    #-- to check
    df_model_nz = df_model %>% dplyr::filter(nonzero==1)
    df_model_z  = df_model %>% dplyr::filter(nonzero==0)
    dim(df_model_nz)[1]/dim(df_model_z)[1] # ~1?
    
    # All parameter names (exactly as in your model)
    param_names = c(
      "th_ROS_microbe",
      "th_ROS_epith_recover",
      "epith_recovery_chance",
      "rat_com_pat_threshold",
      "diffusion_speed_DAMPs",
      "diffusion_speed_SAMPs",
      "diffusion_speed_ROS",
      "add_ROS",
      "add_DAMPs",
      "add_SAMPs",
      "ros_decay",
      "DAMPs_decay",
      "SAMPs_decay",
      "activation_threshold_DAMPs",
      "activation_threshold_SAMPs",
      "activity_engulf_M0_baseline",
      "activity_engulf_M1_baseline",
      "activity_engulf_M2_baseline",
      "activity_ROS_M1_baseline",
      "rate_leak_commensal_injury",
      "rate_leak_pathogen_injury",
      "rate_leak_commensal_baseline",
      "active_age_limit",
      "treg_discrimination_efficiency"
    )
    
    # Output folder for saving plots
    outdir = "/Users/burcutepekule/Desktop/tregs/gam_histograms"
    dir.create(outdir, showWarnings = FALSE)
    
    # Function to create one histogram plot for a given parameter
    plot_param_distribution = function(param) {
      ggplot(df_model, aes(x = .data[[param]], fill = factor(nonzero))) +
        geom_density(alpha = 0.4, color = NA) +
        scale_fill_manual(values = c("0" = "gray70", "1" = "#3182bd"),
                          labels = c("Zero outcome", "Nonzero outcome")) +
        labs(
          x = param,
          y = "Density",
          fill = "Outcome",
          title = paste("Distribution of", param, "by outcome (nonzero)")
        ) +
        theme_minimal(base_size = 13)
    }
    
    # # Loop through all parameters, plot and save
    # walk(param_names, function(param) {
    #   p = plot_param_distribution(param)
    #   ggsave(
    #     filename = file.path(outdir, paste0(param, "_histogram.png")),
    #     plot = p,
    #     width = 6,
    #     height = 4,
    #     dpi = 300
    #   )
    #   message("✅ Saved: ", param)
    # })
    
    table(df_model$nonzero)
    
    df_agg = df_model %>%
      dplyr::group_by(param_set_id, across(all_of(param_names))) %>%
      dplyr::summarise(
        nonzero_n = sum(nonzero == 1),
        zero_n    = sum(nonzero == 0),
        .groups="drop"
      )
    
    df_agg = df_agg %>% dplyr::mutate(sum_n=nonzero_n+zero_n)
    unique(df_agg$sum_n) # all 10
    
    # This naturally handles different numbers of replicates per param set.
    # If there’s extra-binomial variation, try family = quasibinomial as a robustness check.
    
    # models the probability that a given parameter set produces a nonzero outcome across its replicates.
    gam_binom_counts = gam(
      cbind(nonzero_n, zero_n) ~
        s(th_ROS_microbe, bs="cs") +
        s(th_ROS_epith_recover, bs="cs") +
        s(epith_recovery_chance, bs="cs") +
        s(rat_com_pat_threshold, bs="cs") +
        s(diffusion_speed_DAMPs, bs="cs") +
        s(diffusion_speed_SAMPs, bs="cs") +
        s(diffusion_speed_ROS, bs="cs") +
        s(add_ROS, bs="cs") + s(add_DAMPs, bs="cs") + s(add_SAMPs, bs="cs") +
        s(ros_decay, bs="cs") + s(DAMPs_decay, bs="cs") + s(SAMPs_decay, bs="cs") +
        s(activation_threshold_DAMPs, bs="cs") + s(activation_threshold_SAMPs, bs="cs") +
        s(activity_engulf_M0_baseline, bs="cs") + s(activity_engulf_M1_baseline, bs="cs") +
        s(activity_engulf_M2_baseline, bs="cs") + s(activity_ROS_M1_baseline, bs="cs") +
        s(rate_leak_commensal_injury, bs="cs") + s(rate_leak_pathogen_injury, bs="cs") +
        s(rate_leak_commensal_baseline, bs="cs") + s(active_age_limit, bs="cs") + #active_age_limit takes integer values from 1 to 30 and has a natural ordering (larger numbers mean something systematically different, not just a label), it is better treated as numeric — not as a factor.
        s(treg_discrimination_efficiency, bs="cs"),
      data   = df_agg,
      family = binomial,
      method = "REML"
    )
    
    saveRDS(gam_binom_counts, paste0('gam_binom_counts_',ss_start_threshold,"_",inj_type,'.rds'))
    summary(gam_binom_counts)      # see which predictors matter
    plot(gam_binom_counts, pages=1, shade=TRUE)   # visualize smooth effects - log odds!
    gam.check(gam_binom_counts)    # residual diagnostics
    dev.off()
    plot(gam_binom_counts, select = 3, shade = TRUE, seWithMean = TRUE, rug = TRUE)
    # Use gratia package
    library(gratia)
    # dev.off()
    # png(paste0(outdir,"/gam_effects_probs.png"), width = 12, height = 10, units = "in", res = 300)
    # gratia::draw(gam_binom_counts, fun = plogis)
    # dev.off()
    
    # ----------- Systematic Approach for Non-Monotonic Effects: Find Regions!
    
    # When you use fun = plogis in gratia::draw(), it transforms the partial effect
    # (the smooth term alone) to a probability-like scale, but this isn't
    # the same as the actual probability because it doesn't include the intercept
    # and other terms. However, for finding regions, this doesn't matter! We
    # just need to find where the effect is "high" relative to itself.
    
    # Method 1: Extract data from gratia plot object
    extract_regions_from_plot = function(gam_model, param_name, prob_threshold = 0.5) {
      
      # Generate the plot with plogis transformation
      p = draw(gam_model, select = paste0("s(", param_name, ")"), fun = plogis, all.terms = TRUE)
      
      # Extract the data layer
      plot_data = layer_data(p, 2)  # Get first layer (the line)
      
      # The data has x (parameter values) and y (transformed effect)
      result_data = data.frame(
        param_value = plot_data$x,
        effect = plot_data$y
      )
      
      # Find regions above threshold
      result_data = result_data %>%
        arrange(param_value) %>%
        mutate(
          above_threshold = effect >= prob_threshold,
          region_change = above_threshold != lag(above_threshold, default = FALSE),
          region_id = cumsum(region_change)
        )
      
      # Extract region bounds
      regions = result_data %>%
        filter(above_threshold) %>%
        group_by(region_id) %>%
        summarise(
          lower = min(param_value),
          upper = max(param_value),
          mean_effect = mean(effect),
          max_effect = max(effect),
          .groups = "drop"
        )
      
      return(list(
        parameter = param_name,
        n_regions = nrow(regions),
        regions = regions,
        plot_data = result_data
      ))
    }
    
    # Test with th_ROS_microbe
    th_ros_regions = extract_regions_from_plot(gam_binom_counts, "th_ROS_microbe", prob_threshold = 0.5)
    
    cat("th_ROS_microbe has", th_ros_regions$n_regions, "region(s):\n")
    print(th_ros_regions$regions)
    
    # # Verify visually
    # p = draw(gam_binom_counts, select = "s(th_ROS_microbe)", fun = plogis) +
    #   geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")
    # 
    # if (nrow(th_ros_regions$regions) > 0) {
    #   for (i in 1:nrow(th_ros_regions$regions)) {
    #     p = p + annotate("rect",
    #                       xmin = th_ros_regions$regions$lower[i],
    #                       xmax = th_ros_regions$regions$upper[i],
    #                       ymin = -Inf, ymax = Inf,
    #                       alpha = 0.2, fill = "green")
    #   }
    # }
    # print(p)
    
    ### loop over all?
    # Process all parameters
    all_param_regions = lapply(param_names, function(p) {
      tryCatch({
        extract_regions_from_plot(gam_binom_counts, p, prob_threshold = 0.5)
      }, error = function(e) {
        message("Error processing ", p, ": ", e$message)
        list(parameter = p, n_regions = 0, regions = data.frame(), plot_data = NULL)
      })
    })
    names(all_param_regions) = param_names
    
    # Create summary dataframe
    region_summary = data.frame(
      parameter = param_names,
      n_regions = sapply(all_param_regions, function(x) x$n_regions)
    ) %>%
      arrange(desc(n_regions))
    
    print(region_summary)
    
    # Create detailed bounds dataframe
    all_bounds = do.call(rbind, lapply(param_names, function(p) {
      regions = all_param_regions[[p]]$regions
      if (nrow(regions) > 0) {
        regions %>%
          mutate(
            parameter = p,
            region_number = row_number()
          ) %>%
          select(parameter, region_number, lower, upper, mean_effect, max_effect)
      } else {
        data.frame(
          parameter = p,
          region_number = 0,
          lower = NA,
          upper = NA,
          mean_effect = NA,
          max_effect = NA
        )
      }
    }))
    
    print(all_bounds)
    all_bounds$lower=round(all_bounds$lower,2)
    all_bounds$upper=round(all_bounds$upper,2)
    
    # merge with significance
    df_significance = summary(gam_binom_counts)$s.table %>%
      as.data.frame() %>%
      tibble::rownames_to_column("term") %>%
      mutate(
        parameter = gsub("s\\((.*)\\)", "\\1", term)
      ) %>%
      select(parameter, edf, `Ref.df`, `Chi.sq`, `p-value`)
    
    # Merge with all_bounds
    all_bounds = all_bounds %>%
      left_join(df_significance, by = "parameter")
    
    all_bounds$`p-value`=round(all_bounds$`p-value`,3)
    
    significant_bounds = all_bounds %>%
      filter(`p-value` < sign_treshold, !is.na(lower)) %>%  # Significant and has actual regions
      arrange(`p-value`)
    
    # Save the bounds to CSV
    # write.csv(all_bounds, "parameter_constraint_bounds.csv", row.names = FALSE)
    # write.csv(region_summary, "parameter_region_summary.csv", row.names = FALSE)
    
    # Create and save individual plots
    # for (param in param_names) {
    #   
    #   region_info = all_param_regions[[param]]
    #   
    #   # Create the plot
    #   p = draw(gam_binom_counts, select = paste0("s(", param, ")"), fun = plogis) +
    #     geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.8) +
    #     labs(title = paste0(param, " (", region_info$n_regions, " region", 
    #                         ifelse(region_info$n_regions != 1, "s", ""), ")")) +
    #     theme_minimal(base_size = 12)
    #   
    #   sig_info = all_bounds %>% filter(parameter==param)
    #   p_val = unique(sig_info$`p-value`)
    #   
    #   if(p_val<sign_treshold){
    #     # Add shaded regions
    #     if (nrow(region_info$regions) > 0) {
    #       for (i in 1:nrow(region_info$regions)) {
    #         p = p + annotate("rect",
    #                           xmin = region_info$regions$lower[i],
    #                           xmax = region_info$regions$upper[i],
    #                           ymin = -Inf, ymax = Inf,
    #                           alpha = 0.2, fill = "green")
    #       }
    #     }
    #   }else{
    #     # Add shaded regions
    #     if (nrow(region_info$regions) > 0) {
    #       for (i in 1:nrow(region_info$regions)) {
    #         p = p + annotate("rect",
    #                           xmin = region_info$regions$lower[i],
    #                           xmax = region_info$regions$upper[i],
    #                           ymin = -Inf, ymax = Inf,
    #                           alpha = 0.2, fill = "yellow")
    #       }
    #     }
    #   }
    #   
    #   # Save individual plot
    #   ggsave(
    #     filename = paste0(outdir,"/gam_constraint_plots/", param, ".png"),
    #     plot = p,
    #     width = 8,
    #     height = 6,
    #     dpi = 300
    #   )
    #   
    #   cat("Saved plot for", param, "\n")
    # }
    
    # Create a combined plot with all parameters (grid view)
    plot_list = lapply(param_names, function(param) {
      
      region_info = all_param_regions[[param]]
      
      p = draw(gam_binom_counts, select = paste0("s(", param, ")"), fun = plogis) +
        geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.5) +
        labs(title = paste0(param, " (", region_info$n_regions, ")")) +
        theme_minimal(base_size = 8) +
        theme(
          plot.title = element_text(size = 8),
          axis.title = element_text(size = 7),
          axis.text = element_text(size = 6)
        )
      
      
      sig_info = all_bounds %>% filter(parameter==param)
      p_val = unique(sig_info$`p-value`)
      
      if(p_val<sign_treshold){
        # Add shaded regions
        if (nrow(region_info$regions) > 0) {
          for (i in 1:nrow(region_info$regions)) {
            p = p + annotate("rect",
                             xmin = region_info$regions$lower[i],
                             xmax = region_info$regions$upper[i],
                             ymin = -Inf, ymax = Inf,
                             alpha = 0.2, fill = "green")
          }
        }
      }else{
        # Add shaded regions
        if (nrow(region_info$regions) > 0) {
          for (i in 1:nrow(region_info$regions)) {
            p = p + annotate("rect",
                             xmin = region_info$regions$lower[i],
                             xmax = region_info$regions$upper[i],
                             ymin = -Inf, ymax = Inf,
                             alpha = 0.2, fill = "yellow")
          }
        }
      }
      return(p)
    })
    library(gridExtra)
    
    # Save combined plot
    combined_plot = grid.arrange(grobs = plot_list, ncol = 5)
    ggsave(
      filename = paste0(outdir,"/all_parameters_combined_",ss_start_threshold,"_",inj_type,".png"),
      plot = combined_plot,
      width = 10,
      height = 7,
      dpi = 300
    )
    
    # # Loop through all parameters, plot and save
    # walk(param_names, function(param) {
    #   p = plot_param_distribution(param)
    #   new_bounds_var = all_bounds %>% filter(parameter==param)
    #   
    #   for(i in new_bounds_var$region_number){
    #     lo = new_bounds_var[i,]$lower
    #     hi = new_bounds_var[i,]$upper
    #     p  = p +
    #       geom_vline(xintercept = lo, color = "red", linetype = "dashed", linewidth = 0.6) +
    #       geom_vline(xintercept = hi, color = "red", linetype = "dashed", linewidth = 0.6)
    #   }
    # 
    #   ggsave(
    #     filename = file.path(outdir, paste0(param, "_histogram.png")),
    #     plot = p,
    #     width = 6,
    #     height = 4,
    #     dpi = 300
    #   )
    #   message("✅ Saved with bounds: ", param, " [", round(lo,3), ", ", round(hi,3), "]")
    # })
    
    
    # Create a combined plot with all parameter histograms (grid view)
    hist_plot_list = lapply(param_names, function(param) {
      
      # Create base histogram
      p = ggplot(df_model, aes(x = .data[[param]], fill = factor(nonzero))) +
        geom_density(alpha = 0.4, bw = 0.1) +
        scale_fill_manual(values = c("0" = "gray70", "1" = "#3182bd"),
                          labels = c("Zero outcome", "Nonzero outcome")) +
        labs(
          x = param,
          y = "Density",
          fill = "Outcome"
        ) +
        theme_minimal(base_size = 7) +
        theme(
          plot.title = element_text(size = 8),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 5),
          legend.position = "none"  # Remove individual legends for cleaner combined plot
        )
      
      # Add vertical lines for constraint bounds
      new_bounds_var = all_bounds %>% filter(parameter == param)
      
      if (nrow(new_bounds_var) > 0) {
        for (i in new_bounds_var$region_number) {
          lo = new_bounds_var[i, ]$lower
          hi = new_bounds_var[i, ]$upper
          if (!is.na(lo) && !is.na(hi)) {
            p = p +
              geom_vline(xintercept = lo, color = "red", linetype = "dashed", linewidth = 0.4) +
              geom_vline(xintercept = hi, color = "red", linetype = "dashed", linewidth = 0.4)
          }
        }
      }
      
      return(p)
    })
    
    # Create a separate legend plot
    legend_plot = ggplot(df_model, aes(x = th_ROS_microbe, fill = factor(nonzero))) +
      geom_density(alpha = 0.4) +  
      scale_fill_manual(values = c("0" = "gray70", "1" = "#3182bd"),
                        labels = c("Zero outcome", "Nonzero outcome")) +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 9))
    
    # Extract the legend
    library(ggpubr)
    legend = get_legend(legend_plot)
    
    
    # Combine histograms
    combined_hist_plot = grid.arrange(
      grobs = hist_plot_list, 
      ncol = 5,
      bottom = legend
    )
    
    # Save combined histogram plot
    ggsave(
      filename = file.path(outdir, paste0("/all_parameters_histograms_combined_",ss_start_threshold,"_",inj_type,".png")),
      plot = combined_hist_plot,
      width = 10,
      height = 7,
      dpi = 300
    )
  }
}
