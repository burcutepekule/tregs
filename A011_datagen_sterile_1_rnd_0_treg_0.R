rm(list=ls())
library(dplyr)
library(tidyr)
# library(cowplot)
# library(ggplot2)
# library(gridExtra)
# library(grid)
# library(av)

split_equal = function(x, n_chunks) {
  split(x, cut(seq_along(x), breaks = n_chunks, labels = FALSE))
}
args = commandArgs(trailingOnly = TRUE)
n1   = as.integer(args[1])
n2   = as.integer(args[2])

chunks    = split_equal(0:9999, n1)
loop_over = chunks[[n2]]

dir_name = './frames'
dir.create(dir_name, showWarnings = FALSE)

dir_name_data = './mass_sim_results_R'
dir.create(dir_name_data, showWarnings = FALSE)

# random_stream_file = paste0("./random_streams/random_numbers_seed_0.csv")
# stream_in          = scan(random_stream_file, quiet = TRUE, skip = 1)

# Default values
sterile              = 1   # 0 = infection, 1 = sterile injury
allow_tregs          = 0   # Allow tregs to do their job
randomize_tregs      = 0   # 0 = follow DAMPs, 1 = random movement
use_synchronized_rng = TRUE  # Use synchronized random numbers for fair comparisons

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

# ============================================================================
# READ PARAMETERS FROM CSV
# ============================================================================
params_df = read.csv("./original_lhs_parameters.csv", stringsAsFactors = FALSE)

# ============================================================================
# FIXED PARAMETERS (not in CSV)
# ============================================================================
t_max      = 5000
plot_on    = 0
plot_every = Inf
grid_size  = 25
n_phagocytes = round(grid_size * grid_size * 0.35)
n_tregs = round(grid_size * grid_size * 0.35)
n_commensals_lp = 20

injury_percentage = 60
max_level_injury = 5

max_cell_value_ROS = 1
max_cell_value_DAMPs = 1
max_cell_value_SAMPs = 1

lim_ROS = max_cell_value_ROS
lim_DAMP = max_cell_value_DAMPs
lim_SAMP = max_cell_value_SAMPs

act_radius_ROS = 1
act_radius_treg = 1
act_radius_DAMPs = 1
act_radius_SAMPs = 1

# Logistic function parameters (for epithelial injury calculation)
k_in = 0.044
x0_in = 50
shift_by = 10

# for(param_set_id_use in 0:9999){
for(param_set_id_use in loop_over){
  param_set_use = params_df %>% dplyr::filter(param_set_id==param_set_id_use)
  longitudinal_df_keep = c()
  print(paste0('Processing param set ',param_set_id_use,' 😱'))
  
  random_stream_file = paste0("./random_streams/random_numbers_seed_",param_set_id_use,".csv")
  stream_in          = scan(random_stream_file, quiet = TRUE, skip = 1)
  
  for (reps_in in 0:9){
    rm(list = setdiff(ls(),c("param_set_id_use","reps_in","stream_in",
                             "sterile","allow_tregs","randomize_tregs",
                             "use_synchronized_rng","params_df","param_set_use",
                             "colnames_insert","longitudinal_df_keep",
                             "dir_name_data","dir_name","t_max","plot_on",
                             "plot_every","grid_size","n_phagocytes","n_tregs",
                             "n_commensals_lp","injury_percentage","max_level_injury",
                             "max_cell_value_ROS","max_cell_value_DAMPs","max_cell_value_SAMPs",
                             "lim_ROS","lim_DAMP","lim_SAMP","act_radius_ROS","act_radius_treg",
                             "act_radius_DAMPs","act_radius_SAMPs","k_in","x0_in","shift_by",
                             "random_stream_file","stream_in")))
    
    # ============================================================================
    # LOAD FUNCTIONS
    # ============================================================================
    
    source("./MISC/FAST_FUNCTIONS.R")
    source("./MISC/PLOT_FUNCTIONS.R")
    
    # ============================================================================
    # DEAL WITH RNG
    # ============================================================================

    # # shift for every rep, can't afford to have so many rnd generators
    stream_in_use      = c(stream_in[(shift_by*reps_in+1):length(stream_in)], stream_in[1:(shift_by*reps_in)])
    
    # # shift for every rep, can't afford to have so many rnd generators
    # stream_in_use = c(stream_in[(shift_by*reps_in+1):length(stream_in)], stream_in[1:(shift_by*reps_in)])
    
    rng_env        = new.env()
    rng_env$stream = stream_in_use
    rng_env$index  = 1
    
    # Custom uniform sampler
    my_runif = function(n = 1) {
      if (is.null(rng_env$index)) rng_env$index = 1
      
      if (rng_env$index + n - 1 > length(rng_env$stream)) {
        stop(sprintf("Random stream exhausted! Used %d numbers, need %d more.",
                     rng_env$index - 1, n))
      }
      
      vals = rng_env$stream[rng_env$index:(rng_env$index + n - 1)]
      rng_env$index = rng_env$index + n
      vals
    }
    
    # Custom sample function
    my_sample = function(x, size, replace = FALSE, prob = NULL) {
      n = length(x)
      if (!replace && size > n)
        stop("cannot take a sample larger than the population")
      
      probs = if (is.null(prob)) rep(1/n, n) else prob / sum(prob)
      u = my_runif(size)
      cumprob = cumsum(probs)
      indices = findInterval(u, cumprob) + 1
      x[indices]
    }
    
    # Custom beta sampler
    my_rbeta = function(n, shape1, shape2) {
      u1 = my_runif(n)
      u2 = my_runif(n)
      qbeta(u1, shape1, shape2)
    }
    
    # Custom gamma sampler
    my_rgamma = function(n, shape, rate = 1) {
      u = my_runif(n)
      qgamma(u, shape = shape, rate = rate)
    }
    
    # Override base R functions
    runif = my_runif
    sample = my_sample
    rbeta = my_rbeta
    rgamma = my_rgamma
    
    # ============================================================================
    # ASSIGN PARAMETERS FROM CSV
    # ============================================================================
    # Thresholds
    th_ROS_microbe = param_set_use$th_ROS_microbe
    th_ROS_epith_recover = param_set_use$th_ROS_epith_recover
    epith_recovery_chance = param_set_use$epith_recovery_chance
    rat_com_pat_threshold = param_set_use$rat_com_pat_threshold
    
    # Diffusion speeds
    diffusion_speed_DAMPs = param_set_use$diffusion_speed_DAMPs
    diffusion_speed_SAMPs = param_set_use$diffusion_speed_SAMPs
    diffusion_speed_ROS = param_set_use$diffusion_speed_ROS
    
    # Signal production
    add_ROS = param_set_use$add_ROS
    add_DAMPs = param_set_use$add_DAMPs
    add_SAMPs = param_set_use$add_SAMPs
    
    # Decay rates
    ros_decay = param_set_use$ros_decay
    DAMPs_decay = param_set_use$DAMPs_decay
    SAMPs_decay = param_set_use$SAMPs_decay
    
    # Activation thresholds
    activation_threshold_DAMPs = param_set_use$activation_threshold_DAMPs
    activation_threshold_SAMPs = param_set_use$activation_threshold_SAMPs
    
    # Engulfment activities
    activity_engulf_M0_baseline = param_set_use$activity_engulf_M0_baseline
    activity_engulf_M1_baseline = param_set_use$activity_engulf_M1_baseline
    activity_engulf_M2_baseline = param_set_use$activity_engulf_M2_baseline
    activity_engulf_max = 0.99
    
    # ROS production activities
    activity_ROS_M0_baseline = 0.00
    activity_ROS_M1_baseline = param_set_use$activity_ROS_M1_baseline
    activity_ROS_M2_baseline = 0.00
    activity_ROS_max = 0.99
    
    # Leak rates
    rate_leak_commensal_injury = param_set_use$rate_leak_commensal_injury
    rate_leak_pathogen_injury = ifelse(sterile == 1, 0.0, param_set_use$rate_leak_pathogen_injury)
    rate_leak_commensal_baseline = param_set_use$rate_leak_commensal_baseline
    
    # Phagocyte parameters
    active_age_limit = as.integer(param_set_use$active_age_limit)
    cc_phagocyte = 5
    digestion_time = 1
    
    # Treg parameters
    treg_vicinity_effect = 1
    treg_discrimination_efficiency = param_set_use$treg_discrimination_efficiency
    # allow_tregs_to_suppress_cognate = FALSE
    
    # ============================================================================
    # INITIALIZE SIMULATION
    # ============================================================================
    injury_site = get_middle_percent(seq(1, grid_size), injury_percentage)
    n_pathogens_lp = round(rate_leak_pathogen_injury * length(injury_site))
    
    precision = 10 * (exp(5 * treg_discrimination_efficiency))
    
    # Calculate step sizes for activity increases
    activity_engulf_M1_step = (activity_engulf_max - activity_engulf_M1_baseline) / cc_phagocyte
    activity_engulf_M2_step = (activity_engulf_max - activity_engulf_M2_baseline) / cc_phagocyte
    activity_ROS_M1_step = (activity_ROS_max - activity_ROS_M1_baseline) / cc_phagocyte
    
    # Initialize fields
    DAMPs = matrix(0, grid_size, grid_size)
    SAMPs = matrix(0, grid_size, grid_size)
    ROS = matrix(0, grid_size, grid_size)
    
    # Initialize epithelium
    epithelium = data.frame(
      x = seq(1, grid_size, 1),
      y = rep(0, grid_size),
      level_injury = 0,
      id = seq(1, grid_size)
    )
    epithelium[injury_site, ]$level_injury = 1
    
    # Initialize phagocytes
    phagocyte_x = sample(1:grid_size, n_phagocytes, TRUE)
    phagocyte_y = sample(2:grid_size, n_phagocytes, TRUE)
    phagocyte_pathogens_engulfed = rep(0, n_phagocytes)
    phagocyte_commensals_engulfed = rep(0, n_phagocytes)
    phagocyte_num_times_activated = rep(0, n_phagocytes)
    phagocyte_phenotype = rep(0, n_phagocytes)  # 0=M0, 1=M1, 2=M2
    phagocyte_activity_ROS = rep(activity_ROS_M0_baseline, n_phagocytes)
    phagocyte_activity_engulf = rep(activity_engulf_M0_baseline, n_phagocytes)
    phagocyte_active_age = rep(0, n_phagocytes)
    phagocyte_bacteria_registry = matrix(0, nrow = n_phagocytes, ncol = cc_phagocyte)
    
    # Initialize tregs
    treg_x = sample(1:grid_size, n_tregs, TRUE)
    treg_y = sample(2:grid_size, n_tregs, TRUE)
    treg_active_age = rep(0, n_tregs)
    treg_phenotype = rep(0, n_tregs)  # 0=resting, 1=activated
    treg_activity_SAMPs_binary = rep(0, n_tregs)
    
    # Initialize pathogens
    if (n_pathogens_lp == 0) {
      pathogen_coords = matrix(numeric(0), ncol = 3)
      colnames(pathogen_coords) = c("x", "y", "id")
    } else {
      pathogen_coords = matrix(c(
        sample(injury_site, n_pathogens_lp, TRUE),
        rep(1, n_pathogens_lp),
        seq(1, n_pathogens_lp)
      ), ncol = 3)
      colnames(pathogen_coords) = c("x", "y", "id")
    }
    
    # Initialize commensals
    commensal_coords = matrix(c(
      sample(1:grid_size, n_commensals_lp, TRUE),
      sample(1:grid_size, n_commensals_lp, TRUE),
      seq(1, n_commensals_lp)
    ), ncol = 3)
    colnames(commensal_coords) = c("x", "y", "id")
    
    last_id_pathogen = n_pathogens_lp
    last_id_commensal = n_commensals_lp
    
    # Death counters
    pathogens_killed_by_ROS = 0
    pathogens_killed_by_Mac = rep(0, 3)
    commensals_killed_by_ROS = 0
    commensals_killed_by_Mac = rep(0, 3)
    
    # Longitudinal tracking
    epithelium_longitudinal = matrix(0, nrow = t_max, ncol = (max_level_injury + 1))
    macrophages_longitudinal = matrix(0, nrow = t_max, ncol = 1 + 2 * (cc_phagocyte + 1))
    microbes_longitudinal = matrix(0, nrow = t_max, ncol = 2)
    tregs_longitudinal = matrix(0, nrow = t_max, ncol = 2)
    microbes_cumdeath_longitudinal = matrix(0, nrow = t_max, ncol = 2 * 4)
    
    # ============================================================================
    # MAIN SIMULATION LOOP
    # ============================================================================
    for (t in 1:t_max) {
      
      # Update injury site
      injury_site_updated = which(epithelium$level_injury > 0)
      
      # ========================================================================
      # UPDATE SAMPs (from activated Tregs)
      # ========================================================================
      active_tregs = which(treg_phenotype == 1)
      if (length(active_tregs) > 0) {
        for (i in active_tregs) { # add *allow_tregs to the end to truly check same randomness
          SAMPs[treg_y[i], treg_x[i]] = SAMPs[treg_y[i], treg_x[i]] +
            treg_activity_SAMPs_binary[i]*add_SAMPs*allow_tregs
        }
      }
      
      # ========================================================================
      # UPDATE ROS (from M1 phagocytes)
      # ========================================================================
      M1_phagocytes = which(phagocyte_phenotype == 1)
      if (length(M1_phagocytes) > 0) {
        for (i in M1_phagocytes) {
          ROS[phagocyte_y[i], phagocyte_x[i]] = ROS[phagocyte_y[i], phagocyte_x[i]] +
            phagocyte_activity_ROS[i] * add_ROS
        }
      }
      
      # ========================================================================
      # MOVE MICROBES (UPDATED ORDER: before DAMP calculation)
      # ========================================================================
      if (nrow(pathogen_coords) > 0) {
        dy = ifelse(pathogen_coords[, "y"] == 1,
                    sample(c(1), size = nrow(pathogen_coords), replace = TRUE),
                    sample(c(-1, 0, 1), size = nrow(pathogen_coords), replace = TRUE))
        dx = iszero_coordinates(dy)
        pathogen_coords[, "x"] = pmin(pmax(pathogen_coords[, "x"] + dx, 1), grid_size)
        pathogen_coords[, "y"] = pmin(pmax(pathogen_coords[, "y"] + dy, 1), grid_size)
      }
      
      if (nrow(commensal_coords) > 0) {
        dy = ifelse(commensal_coords[, "y"] == 1,
                    sample(c(1), size = nrow(commensal_coords), replace = TRUE),
                    sample(c(-1, 0, 1), size = nrow(commensal_coords), replace = TRUE))
        dx = iszero_coordinates(dy)
        commensal_coords[, "x"] = pmin(pmax(commensal_coords[, "x"] + dx, 1), grid_size)
        commensal_coords[, "y"] = pmin(pmax(commensal_coords[, "y"] + dy, 1), grid_size)
      }
      
      # ========================================================================
      # PRE-CALCULATE MICROBE COUNTS AT EPITHELIUM (y=1)
      # ========================================================================
      pathogen_epithelium_counts = rep(0, grid_size)
      if (nrow(pathogen_coords) > 0) {
        epithelium_pathogens = pathogen_coords[pathogen_coords[, "y"] == 1, , drop = FALSE]
        if (nrow(epithelium_pathogens) > 0) {
          pathogen_epithelium_counts = tabulate(epithelium_pathogens[, "x"], nbins = grid_size)
        }
      }
      
      commensal_epithelium_counts = rep(0, grid_size)
      if (nrow(commensal_coords) > 0) {
        epithelium_commensals = commensal_coords[commensal_coords[, "y"] == 1, , drop = FALSE]
        if (nrow(epithelium_commensals) > 0) {
          commensal_epithelium_counts = tabulate(epithelium_commensals[, "x"], nbins = grid_size)
        }
      }
      
      # ========================================================================
      # UPDATE DAMPs (UPDATED: includes epithelium + microbe counts)
      # ========================================================================
      # DAMPs from epithelial injury
      DAMPs[1, ] = DAMPs[1, ] + epithelium$level_injury * add_DAMPs
      
      # DAMPs from microbes touching epithelium (NEW from STANDALONE)
      DAMPs[1, ] = DAMPs[1, ] + 1 * logistic_scaled_0_to_5_quantized(
        pathogen_epithelium_counts + commensal_epithelium_counts
      ) * add_DAMPs
      
      # Additional DAMPs from pathogen locations (currently multiplied by 0 in STANDALONE)
      DAMPs_add = matrix(0, nrow = nrow(DAMPs), ncol = ncol(DAMPs))
      if (nrow(pathogen_coords) > 0) {
        pat_counts = table(pathogen_coords[, "x"], pathogen_coords[, "y"])
        if (dim(pat_counts)[1] > 0) {
          pat_df = as.data.frame(pat_counts)
          names(pat_df) = c("x", "y", "count")
          pat_df$x = as.numeric(as.character(pat_df$x))
          pat_df$y = as.numeric(as.character(pat_df$y))
          pat_df$val = 0 * add_DAMPs * logistic_scaled_0_to_5_quantized(pat_df$count)
          DAMPs_add[cbind(pat_df$y, pat_df$x)] = pat_df$val
        }
      }
      DAMPs = DAMPs + DAMPs_add
      
      # ========================================================================
      # DIFFUSE & DECAY SIGNALS
      # ========================================================================
      DAMPs = diffuse_matrix(DAMPs, diffusion_speed_DAMPs, max_cell_value_DAMPs)
      SAMPs = diffuse_matrix(SAMPs, diffusion_speed_SAMPs, max_cell_value_SAMPs)
      ROS = diffuse_matrix(ROS, diffusion_speed_ROS, max_cell_value_ROS)
      DAMPs = DAMPs - DAMPs_decay * DAMPs
      SAMPs = SAMPs - SAMPs_decay * SAMPs
      ROS = ROS - ros_decay * ROS
      
      # ========================================================================
      # MOVE PHAGOCYTES AND TREGS
      # ========================================================================
      density_matrix_tregs = if (randomize_tregs == 1) 0 * DAMPs else DAMPs
      density_matrix_phagocytes = DAMPs
      
      all_equal_treg = all(density_matrix_tregs == density_matrix_tregs[1, 1])
      all_equal_phagocytes = all(density_matrix_phagocytes == density_matrix_phagocytes[1, 1])
      
      # Move tregs
      if (!all_equal_treg) {
        for (i in 1:length(treg_x)) {
          x = treg_x[i]
          y = treg_y[i]
          x_range = max(1, x - 1):min(grid_size, x + 1)
          y_range = max(1, y - 1):min(grid_size, y + 1)
          neighbors_x = rep(x_range, each = length(y_range))
          neighbors_y = rep(y_range, times = length(x_range))
          neighbor_densities = density_matrix_tregs[cbind(neighbors_y, neighbors_x)]
          total = sum(neighbor_densities)
          if (total > 0) {
            probs = neighbor_densities / total
          } else {
            probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
          }
          chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
          treg_x[i] = neighbors_x[chosen_idx]
          treg_y[i] = neighbors_y[chosen_idx]
        }
      } else {
        dy_treg = ifelse(treg_y == 1,
                         sample(c(1), size = length(treg_y), replace = TRUE),
                         sample(c(-1, 0, 1), size = length(treg_y), replace = TRUE))
        dx_treg = iszero_coordinates(dy_treg)
        treg_x = pmin(pmax(treg_x + dx_treg, 1), grid_size)
        treg_y = pmin(pmax(treg_y + dy_treg, 1), grid_size)
      }
      
      # Move phagocytes
      if (!all_equal_phagocytes) {
        for (i in 1:length(phagocyte_x)) {
          x = phagocyte_x[i]
          y = phagocyte_y[i]
          x_range = max(1, x - 1):min(grid_size, x + 1)
          y_range = max(1, y - 1):min(grid_size, y + 1)
          neighbors_x = rep(x_range, each = length(y_range))
          neighbors_y = rep(y_range, times = length(x_range))
          neighbor_densities = density_matrix_phagocytes[cbind(neighbors_y, neighbors_x)]
          total = sum(neighbor_densities)
          if (total > 0) {
            probs = neighbor_densities / total
          } else {
            probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
          }
          chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
          phagocyte_x[i] = neighbors_x[chosen_idx]
          phagocyte_y[i] = neighbors_y[chosen_idx]
        }
      } else {
        dy_phagocyte = ifelse(phagocyte_y == 1,
                              sample(c(1), size = length(phagocyte_y), replace = TRUE),
                              sample(c(-1, 0, 1), size = length(phagocyte_y), replace = TRUE))
        dx_phagocyte = iszero_coordinates(dy_phagocyte)
        phagocyte_x = pmin(pmax(phagocyte_x + dx_phagocyte, 1), grid_size)
        phagocyte_y = pmin(pmax(phagocyte_y + dy_phagocyte, 1), grid_size)
      }
      
      # ========================================================================
      # ADD NEW MICROBES
      # ========================================================================
      n_pathogens_lp_new = round(mean(epithelium$level_injury) * rate_leak_pathogen_injury *
                                   length(injury_site_updated))
      if (n_pathogens_lp_new > 0) {
        new_pathogen_coords = matrix(c(
          sample(1:grid_size, n_pathogens_lp_new, replace = TRUE, prob = epithelium$level_injury),
          rep(1, n_pathogens_lp_new),
          last_id_pathogen + seq(1, n_pathogens_lp_new)
        ), ncol = 3)
        colnames(new_pathogen_coords) = c("x", "y", "id")
        pathogen_coords = rbind(pathogen_coords, new_pathogen_coords)
        last_id_pathogen = last_id_pathogen + n_pathogens_lp_new
      }
      
      n_commensals_lp_new_injury = round(mean(epithelium$level_injury) * rate_leak_commensal_injury *
                                           length(injury_site_updated))
      n_commensals_lp_new_baseline = round(rate_leak_commensal_baseline * grid_size)
      
      total_new_commensals = n_commensals_lp_new_baseline + n_commensals_lp_new_injury
      if (total_new_commensals > 0) {
        baseline_x = sample(1:grid_size, n_commensals_lp_new_baseline, TRUE)
        injury_x = if (n_commensals_lp_new_injury > 0) {
          sample(1:grid_size, n_commensals_lp_new_injury, TRUE, prob = epithelium$level_injury)
        } else {
          numeric(0)
        }
        new_commensal_coords = matrix(c(
          c(baseline_x, injury_x),
          rep(1, total_new_commensals),
          last_id_commensal + seq(1, total_new_commensals)
        ), ncol = 3)
        colnames(new_commensal_coords) = c("x", "y", "id")
        commensal_coords = rbind(commensal_coords, new_commensal_coords)
        last_id_commensal = last_id_commensal + total_new_commensals
      }
      
      # ========================================================================
      # PLOTTING
      # ========================================================================
      if (plot_on == 1 & (t %% plot_every == 0 | t == 1)) {
        source('./MISC/CONVERT_TO_DATAFRAME.R')
        p = plot_simtime_simple()
        ggsave(
          paste0(dir_name, "/frame_param_",param_set_id_use,"_rep_", reps_in, "_STERILE_", sterile, "_TREGS_",
                 allow_tregs, "_trnd_", randomize_tregs, "_", t, ".png"),
          plot = p,
          width = 12,
          height = 10,
          dpi = 600,
          bg = "white"
        )
      }
      
      # ========================================================================
      # UPDATE PHAGOCYTE PHENOTYPES
      # ========================================================================
      M0_indices = which(phagocyte_phenotype == 0)
      M1_indices = which(phagocyte_phenotype == 1)
      M2_indices = which(phagocyte_phenotype == 2)
      
      # Registry shifting every digestion_time steps
      if (t %% digestion_time == 0) {
        phagocyte_bacteria_registry = cbind(
          matrix(0, nrow = nrow(phagocyte_bacteria_registry), ncol = 1),
          phagocyte_bacteria_registry[, -ncol(phagocyte_bacteria_registry)]
        )
      }
      
      # Process M0 phagocytes (candidates for activation)
      if (length(M0_indices) > 0) {
        for (i in M0_indices) {
          avg_DAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_DAMPs, DAMPs)
          avg_SAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_SAMPs, SAMPs)
          bacteria_count = sum(phagocyte_bacteria_registry[i, ])
          
          if (avg_DAMPs >= activation_threshold_DAMPs && avg_DAMPs > avg_SAMPs) {
            phagocyte_phenotype[i] = 1
            phagocyte_active_age[i] = 1
            phagocyte_activity_ROS[i] = activity_ROS_M1_baseline + activity_ROS_M1_step * bacteria_count
            phagocyte_activity_engulf[i] = activity_engulf_M1_baseline + activity_engulf_M1_step * bacteria_count
          } else if (avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > avg_DAMPs) {
            phagocyte_phenotype[i] = 2
            phagocyte_active_age[i] = 1
            phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M2_baseline + activity_engulf_M2_step * bacteria_count
          }
        }
      }
      
      # Process M1/M2 phagocytes
      active_indices = c(M1_indices, M2_indices)
      if (length(active_indices) > 0) {
        phagocyte_active_age[active_indices] = phagocyte_active_age[active_indices] + 1
        old_enough = phagocyte_active_age[active_indices] >= active_age_limit
        candidates = active_indices[old_enough]
        
        for (i in candidates) {
          avg_DAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_DAMPs, DAMPs)
          avg_SAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_SAMPs, SAMPs)
          bacteria_count = sum(phagocyte_bacteria_registry[i, ])
          
          if (avg_DAMPs >= activation_threshold_DAMPs && avg_DAMPs > avg_SAMPs) {
            phagocyte_phenotype[i] = 1
            phagocyte_active_age[i] = 1
            phagocyte_activity_ROS[i] = activity_ROS_M1_baseline + activity_ROS_M1_step * bacteria_count
            phagocyte_activity_engulf[i] = activity_engulf_M1_baseline + activity_engulf_M1_step * bacteria_count
          } else if (avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > avg_DAMPs) {
            phagocyte_phenotype[i] = 2
            phagocyte_active_age[i] = 1
            phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M2_baseline + activity_engulf_M2_step * bacteria_count
          } else if (avg_SAMPs < activation_threshold_SAMPs && avg_DAMPs < activation_threshold_DAMPs) {
            phagocyte_phenotype[i] = 0
            phagocyte_active_age[i] = 0
            phagocyte_activity_ROS[i] = activity_ROS_M0_baseline
            phagocyte_activity_engulf[i] = activity_engulf_M0_baseline
          }
        }
      }
      
      # ========================================================================
      # UPDATE TREG ACTIVE AGE
      # ========================================================================
      active_treg_indices = which(treg_phenotype == 1)
      if (length(active_treg_indices) > 0) {
        old_tregs = active_treg_indices[treg_active_age[active_treg_indices] >= active_age_limit]
        young_tregs = active_treg_indices[treg_active_age[active_treg_indices] < active_age_limit]
        
        if (length(young_tregs) > 0) {
          treg_active_age[young_tregs] = treg_active_age[young_tregs] + 1
        }
        
        if (length(old_tregs) > 0) {
          treg_phenotype[old_tregs] = 0
          treg_active_age[old_tregs] = 0
          treg_activity_SAMPs_binary[old_tregs] = 0
        }
      }
      
      # ========================================================================
      # ENGULFMENT PROCESS
      # ========================================================================
      phagocyte_positions = paste(phagocyte_x, phagocyte_y, sep = "_")
      
      for (i in 1:length(phagocyte_x)) {
        px = phagocyte_x[i]
        py = phagocyte_y[i]
        
        # Pathogen engulfment
        if (nrow(pathogen_coords) > 0) {
          pathogen_overlap = (pathogen_coords[, "x"] == px) & (pathogen_coords[, "y"] == py)
          pathogen_indices = which(pathogen_overlap)
          
          if (length(pathogen_indices) > 0) {
            engulf_success = runif(length(pathogen_indices)) < phagocyte_activity_engulf[i]
            indices_to_engulf = pathogen_indices[engulf_success]
            
            if (length(indices_to_engulf) > 0) {
              phagocyte_pathogens_engulfed[i] = phagocyte_pathogens_engulfed[i] + length(indices_to_engulf)
              pathogen_coords = pathogen_coords[-indices_to_engulf, , drop = FALSE]
              phagocyte_bacteria_registry[i, ] = shift_insert_fast(
                phagocyte_bacteria_registry[i, ],
                rep(1, length(indices_to_engulf))
              )
              phagocyte_phenotype_index = phagocyte_phenotype[i] + 1
              pathogens_killed_by_Mac[phagocyte_phenotype_index] =
                pathogens_killed_by_Mac[phagocyte_phenotype_index] + length(indices_to_engulf)
            }
          }
        }
        
        # Commensal engulfment
        if (nrow(commensal_coords) > 0) {
          commensal_overlap = (commensal_coords[, "x"] == px) & (commensal_coords[, "y"] == py)
          commensal_indices = which(commensal_overlap)
          
          if (length(commensal_indices) > 0) {
            engulf_success = runif(length(commensal_indices)) < phagocyte_activity_engulf[i]
            indices_to_engulf = commensal_indices[engulf_success]
            
            if (length(indices_to_engulf) > 0) {
              phagocyte_commensals_engulfed[i] = phagocyte_commensals_engulfed[i] + length(indices_to_engulf)
              commensal_coords = commensal_coords[-indices_to_engulf, , drop = FALSE]
              phagocyte_bacteria_registry[i, ] = shift_insert_fast(
                phagocyte_bacteria_registry[i, ],
                rep(1, length(indices_to_engulf))
              )
              phagocyte_phenotype_index = phagocyte_phenotype[i] + 1
              commensals_killed_by_Mac[phagocyte_phenotype_index] =
                commensals_killed_by_Mac[phagocyte_phenotype_index] + length(indices_to_engulf)
            }
          }
        }
      }
      
      # ========================================================================
      # TREG ACTIVATION & EFFECTOR ACTIONS (UPDATED: Beta distribution sampling)
      # CRITICAL: Always consume random numbers to maintain stream synchronization!
      # ========================================================================
      M1_phagocyte_indices = which(phagocyte_phenotype == 1)
      M2_phagocyte_indices = which(phagocyte_phenotype == 2)
      M_activate_phagocyte_indices = c(M1_phagocyte_indices,M2_phagocyte_indices)
      
      # OLD - ONLY M1 - THINK WHETHER THIS MADE MORE SENSE?
      # if (length(M1_phagocyte_indices) > 0) {
      #   for (i in M1_phagocyte_indices) {
      
      if (length(M_activate_phagocyte_indices) > 0) {
        for (i in M_activate_phagocyte_indices) {
          px = phagocyte_x[i]
          py = phagocyte_y[i]
          
          treg_distances_x = abs(treg_x - px)
          treg_distances_y = abs(treg_y - py)
          nearby_treg_indices = which(treg_distances_x <= treg_vicinity_effect &
                                        treg_distances_y <= treg_vicinity_effect)
          
          if (length(nearby_treg_indices) > 0) {
            num_pat_antigens = phagocyte_pathogens_engulfed[i]
            num_com_antigens = phagocyte_commensals_engulfed[i]
            
            if ((num_pat_antigens + num_com_antigens) > 0) {
              # ALWAYS calculate and sample (to maintain stream synchronization)
              rat_com_pat_real = num_com_antigens / (num_com_antigens + num_pat_antigens)
              alpha = (1 - treg_discrimination_efficiency) * 1 +
                treg_discrimination_efficiency * (rat_com_pat_real * precision)
              beta = (1 - treg_discrimination_efficiency) * 1 +
                treg_discrimination_efficiency * ((1 - rat_com_pat_real) * precision)
              
              # ALWAYS consume from random stream (both Tregs ON and OFF)
              rat_com_pat = sample_rbeta(alpha, beta)
              
              # But ONLY apply the effect if Tregs are allowed to work
              if (rat_com_pat > rat_com_pat_threshold) {
                treg_phenotype[nearby_treg_indices] = 1
                treg_activity_SAMPs_binary[nearby_treg_indices] = 1
                treg_active_age[nearby_treg_indices] = 1
                
                # if (allow_tregs_to_suppress_cognate) {
                #   phagocyte_phenotype[i] = 2
                #   phagocyte_active_age[i] = 1
                #   bacteria_count = sum(phagocyte_bacteria_registry[i, ])
                #   phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
                #   phagocyte_activity_engulf[i] = activity_engulf_M2_baseline +
                #     activity_engulf_M2_step * bacteria_count
                # }
              }
            }
          }
        }
      }
      
      # ========================================================================
      # KILL MICROBES WITH ROS
      # ========================================================================
      if (nrow(pathogen_coords) > 0) {
        pathogen_avg_ROS = numeric(nrow(pathogen_coords))
        for (i in 1:nrow(pathogen_coords)) {
          pathogen_avg_ROS[i] = get_8n_avg_signal_fast(
            pathogen_coords[i, "x"],
            pathogen_coords[i, "y"],
            act_radius_ROS, ROS
          )
        }
        pathogens_to_kill = which(pathogen_avg_ROS > th_ROS_microbe)
        if (length(pathogens_to_kill) > 0) {
          pathogen_coords = pathogen_coords[-pathogens_to_kill, , drop = FALSE]
          pathogens_killed_by_ROS = pathogens_killed_by_ROS + length(pathogens_to_kill)
        }
      }
      
      if (nrow(commensal_coords) > 0) {
        commensal_avg_ROS = numeric(nrow(commensal_coords))
        for (i in 1:nrow(commensal_coords)) {
          commensal_avg_ROS[i] = get_8n_avg_signal_fast(
            commensal_coords[i, "x"],
            commensal_coords[i, "y"],
            act_radius_ROS, ROS
          )
        }
        commensals_to_kill = which(commensal_avg_ROS > th_ROS_microbe)
        if (length(commensals_to_kill) > 0) {
          commensal_coords = commensal_coords[-commensals_to_kill, , drop = FALSE]
          commensals_killed_by_ROS = commensals_killed_by_ROS + length(commensals_to_kill)
        }
      }
      
      # ========================================================================
      # UPDATE EPITHELIAL INJURY (UPDATED: only pathogens, logistic function)
      # ========================================================================
      for (i in 1:nrow(epithelium)) {
        px = epithelium$x[i]
        
        # Get ROS values in vicinity
        x_coordinates = pmax(1, pmin(grid_size, (px - act_radius_ROS):(px + act_radius_ROS)))
        ros_values = ROS[1, x_coordinates]
        mean_ros = mean(ros_values)
        
        # UPDATED: Increase level_injury based on pathogen count ONLY (using logistic function)
        count_pathogens = pathogen_epithelium_counts[px]
        epithelium$level_injury[i] = epithelium$level_injury[i] +
          logistic_scaled_0_to_5_quantized(count_pathogens)
        
        # NOTE: Commensals NO LONGER contribute to epithelial injury (removed from STANDALONE)
        
        # Increase level_injury based on ROS
        if (mean_ros > th_ROS_epith_recover) {
          epithelium$level_injury[i] = epithelium$level_injury[i] + 1
        }
        
        # Apply maximum injury constraint
        epithelium$level_injury[i] = min(epithelium$level_injury[i], max_level_injury)
        
        # RECOVERY: Stochastic recovery when injured
        if (epithelium$level_injury[i] > 0 && runif(1) < epith_recovery_chance) {
          epithelium$level_injury[i] = max(0, epithelium$level_injury[i] - 1)
        }
      }
      
      # ========================================================================
      # SAVE ABUNDANCES
      # ========================================================================
      epithelium_longitudinal[t, ] = as.numeric(table(factor(epithelium$level_injury, levels = 0:5)))
      
      phagocyte_counts = c(
        sum(phagocyte_phenotype == 0),
        tabulate(phagocyte_active_age[phagocyte_phenotype == 1] + 1, cc_phagocyte + 1),
        tabulate(phagocyte_active_age[phagocyte_phenotype == 2] + 1, cc_phagocyte + 1)
      )
      macrophages_longitudinal[t, ] = phagocyte_counts
      
      microbes_longitudinal[t, ] = c(nrow(commensal_coords), nrow(pathogen_coords))
      tregs_longitudinal[t, ] = c(sum(treg_phenotype == 0), sum(treg_phenotype == 1))
      microbes_cumdeath_longitudinal[t, ] = c(
        commensals_killed_by_ROS, commensals_killed_by_Mac,
        pathogens_killed_by_ROS, pathogens_killed_by_Mac
      )
    }
    
    # ============================================================================
    # SAVE LONGITUDINAL DATA
    # ============================================================================
    longitudinal_df = data.frame(
      epithelium_longitudinal,
      macrophages_longitudinal,
      microbes_longitudinal,
      tregs_longitudinal,
      microbes_cumdeath_longitudinal
    )
    
    longitudinal_df$t = 1:t_max
    longitudinal_df$sterile = sterile
    longitudinal_df$tregs_on = allow_tregs
    # longitudinal_df$allow_tregs_to_suppress_cognate = allow_tregs_to_suppress_cognate
    longitudinal_df$randomize_tregs = randomize_tregs
    longitudinal_df$param_set_id = param_set_use$param_set_id
    longitudinal_df$rep_id = reps_in
    
    longitudinal_df = longitudinal_df %>%
      select(t, sterile, tregs_on, #allow_tregs_to_suppress_cognate,
             randomize_tregs, param_set_id, rep_id, everything())
    longitudinal_df_keep = rbind(longitudinal_df_keep, longitudinal_df)
  }
  colnames(longitudinal_df_keep)[c(7:37)] = colnames_insert
  saveRDS(longitudinal_df_keep, paste0(dir_name_data,'/longitudinal_df_param_set_id_',param_set_id_use,
                                       '_sterile_',sterile,
                                       '_trnd_',randomize_tregs,
                                       '_tregs_',allow_tregs,'.rds'))
  print(paste0('Data for param set id ',param_set_id_use,' saved 🥳.'))
}




