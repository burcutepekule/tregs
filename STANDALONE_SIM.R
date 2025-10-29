rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)

# set.seed(42)

rng_env <- new.env()
rng_env$stream <- scan("/Users/burcutepekule/Dropbox/Treg_problem_v2/python/random_numbers_seed_42.txt")
rng_env$index <- 1

# Define the custom uniform sampler
my_runif <- function(n = 1) {
  # check if we have enough numbers left
  if (is.null(rng_env$index)) rng_env$index <- 1
  
  if (rng_env$index + n - 1 > length(rng_env$stream))
    stop("Stream exhausted!")
  
  # extract values
  vals <- rng_env$stream[rng_env$index:(rng_env$index + n - 1)]
  
  # update the pointer
  rng_env$index <- rng_env$index + n
  
  vals
}

my_sample <- function(x, size, replace = FALSE, prob = NULL) {
  n <- length(x)
  if (!replace && size > n)
    stop("cannot take a sample larger than the population")
  
  probs <- if (is.null(prob)) rep(1/n, n) else prob / sum(prob)
  u <- my_runif(size)
  cumprob <- cumsum(probs)
  indices <- findInterval(u, cumprob) + 1
  x[indices]
}

my_rbeta <- function(n, shape1, shape2) {
  u1 <- my_runif(n)
  u2 <- my_runif(n)
  qbeta(u1, shape1, shape2)
}

my_rgamma <- function(n, shape, rate = 1) {
  u <- my_runif(n)
  qgamma(u, shape = shape, rate = rate)
}

# If you really want runif() itself to use your imported stream (so existing code keeps working), 
# you can temporarily mask it in your workspac
runif <- my_runif
sample <- my_sample
rbeta  <- my_rbeta
rgamma = my_rgamma

#---------- SIMULATION FUNCTIONS ------------------------------------------------

sample_rbeta = function(alpha, beta) {
  # Sample from Gamma distributions
  x <- rgamma(1, shape = alpha, rate = 1.0)
  y <- rgamma(1, shape = beta, rate = 1.0)
  
  # Beta is the ratio
  return(x / (x + y))
}
get_middle_percent <- function(seq_vector, percent) {
  n_total <- length(seq_vector)
  n_select <- ceiling(n_total * percent / 100)
  
  # Calculate start and end index for middle values
  mid <- floor(n_total / 2)
  half_window <- floor(n_select / 2)
  
  start_idx <- max(1, mid - half_window + 1)
  end_idx <- min(n_total, start_idx + n_select - 1)
  
  return(seq_vector[start_idx:end_idx])
}

# Optimized version of iszero_coordinates
iszero_coordinates <- function(x) {
  # Initialize with default sampling: -1, 0, 1
  y <- sample(c(-1, 0, 1), length(x), replace = TRUE)
  # Replace where x == 0 with sample from -1 or 1
  zero_idx <- which(x == 0)
  y[zero_idx] <- sample(c(-1, 1), length(zero_idx), replace = TRUE)
  
  return(y)
}

logistic_scaled_0_to_5_quantized <- function(x,  k = k_in, x0 = x0_in) {
  return(round(5*plogis(x, location = x0, scale = 1 / k)))
}

diffuse_matrix <- function(mat, D, max_cell_value) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  # Pad the matrix with zeros around the edges
  padded <- matrix(0, nrow = nr + 2, ncol = nc + 2)
  padded[2:(nr + 1), 2:(nc + 1)] <- mat
  
  # # Compute 4-neighbor diffusion
  # laplacian <- padded[1:nr,   2:(nc+1)] +  # up
  #   padded[3:(nr+2), 2:(nc+1)] +  # down
  #   padded[2:(nr+1), 1:nc] +      # left
  #   padded[2:(nr+1), 3:(nc+2)] -  # right
  #   4 * mat                      # center
  
  # Compute 8-neighbor Laplacian (Moore neighborhood)
  laplacian <- (
    padded[1:nr,     1:nc    ] +  # top-left
      padded[1:nr,     2:(nc+1)] +  # top
      padded[1:nr,     3:(nc+2)] +  # top-right
      padded[2:(nr+1), 1:nc    ] +  # left
      padded[2:(nr+1), 3:(nc+2)] +  # right
      padded[3:(nr+2), 1:nc    ] +  # bottom-left
      padded[3:(nr+2), 2:(nc+1)] +  # bottom
      padded[3:(nr+2), 3:(nc+2)]    # bottom-right
    - 8 * mat                   # center subtraction
  )
  
  mat_new <- mat + D * laplacian
  mat_new <- matrix(pmin(max_cell_value, mat_new), nrow = nrow(mat), ncol = ncol(mat))
  
  return(mat_new)
}

get_8n_avg_signal_fast <- function(x, y, act_radius_signal, signal_matrix) {
  loc  = c(x, y)
  x_coordinates = (loc[1]-act_radius_signal):(loc[1]+act_radius_signal)
  x_coordinates = x_coordinates[x_coordinates>0 & x_coordinates<=grid_size]
  y_coordinates = (loc[2]-act_radius_signal):(loc[2]+act_radius_signal)
  y_coordinates = y_coordinates[y_coordinates>0 & y_coordinates<=grid_size]
  dval = signal_matrix[y_coordinates, x_coordinates]
  return(mean(dval))
}

# Fast shift_insert for matrix operations
shift_insert_fast <- function(vec, insert_vals) {
  n_insert <- length(insert_vals)
  n_vec <- length(vec)
  
  if(n_insert >= n_vec) {
    return(insert_vals[1:n_vec])
  } else {
    return(c(insert_vals, vec[1:(n_vec - n_insert)]))
  }
}
#---------- SIMULATION FUNCTIONS ------------------------------------------------

# Grid parameters
grid_size     = 25
t_max         = 500

# Agent counts 
n_phagocytes      = round(grid_size*grid_size*0.35)
n_tregs           = round(grid_size*grid_size*0.35)

# Initial microbe counts
n_commensals_lp   = 20

# Injury parameters
injury_percentage = 60
max_level_injury  = 5

# Diffusion parameters (must be < 0.125 for stability)
diffusion_speed_DAMPs = 0.07209796
diffusion_speed_SAMPs = 0.08999438
diffusion_speed_ROS   = 0.05232508

# Decay rates
DAMPs_decay = 0.6555613
SAMPs_decay = 0.2551983
ros_decay   = 0.02983629

# Signal production rates
add_DAMPs = 0.4806108
add_SAMPs = 0.8993819
add_ROS   = 0.8972294

# Activation thresholds
activation_threshold_DAMPs = 0.03546454
activation_threshold_SAMPs = 0.7854582

# Killing thresholds
th_ROS_microbe        = 0.2691068
th_ROS_epith_recover  = 0.8875888
epith_recovery_chance = 0.473025

# Engulfment activities (baseline rates)
activity_engulf_M0_baseline = 0.07882888
activity_engulf_M1_baseline = 0.1559077
activity_engulf_M2_baseline = 0.2466745
activity_engulf_max         = 0.99

# ROS production activities
activity_ROS_M0_baseline  = 0.00 
activity_ROS_M1_baseline  = 0.4768265
activity_ROS_M2_baseline  = 0.00 
activity_ROS_max          = 0.99

# Treg parameters
treg_vicinity_effect           = 1 # if 0, that means has to be at the very same pixel
treg_discrimination_efficiency = 0.03508874
rat_com_pat_threshold          = 0.7703095
allow_tregs_to_do_their_job    = 1
allow_tregs_to_suppress_cognate= FALSE
randomize_tregs                = 1

# Phagocyte parameters
cc_phagocyte           = 5
active_age_limit       = 15
digestion_time         = 1

# Infection parameters 
rate_leak_pathogen_injury    = 0.7240775
rate_leak_commensal_injury   = 0.9778595
rate_leak_commensal_baseline = 0.2162822
sterile                      = 0

# Logistic function parameters
k_in  = 0.044
x0_in = 50

# Action radii
act_radius_ROS          = 1 # ROS RADIUS, SMALLER THE LOCAL
act_radius_DAMPs        = 1 # DAMP RADIUS, SMALLER THE LOCAL
act_radius_SAMPs        = 1 # SAMP RADIUS, SMALLER THE LOCAL

injury_site = get_middle_percent(seq(1,grid_size), injury_percentage)
precision   = 10*(exp(5*treg_discrimination_efficiency)) # for later beta sampling

max_cell_value_ROS= 1
max_cell_value_DAMPs = 1
max_cell_value_SAMPs = 1

# for plotting purposes
lim_ROS = max_cell_value_ROS
lim_DAMP= max_cell_value_DAMPs
lim_SAMP= max_cell_value_SAMPs

#-------------------------------ABM SIMULATION ------------------------------------------------------------------------
if(sterile==1){
  rate_leak_pathogen_injury      = 0.00
}else{
  rate_leak_pathogen_injury     = rate_leak_pathogen_injury
}

n_pathogens_lp  = round(rate_leak_pathogen_injury*length(injury_site)) 

# Initialize fields
DAMPs  = matrix(0, grid_size, grid_size)
SAMPs  = matrix(0, grid_size, grid_size)
ROS    = matrix(0, grid_size, grid_size)

activity_engulf_M1_step   = (activity_engulf_max-activity_engulf_M1_baseline)/cc_phagocyte
activity_engulf_M2_step   = (activity_engulf_max-activity_engulf_M2_baseline)/cc_phagocyte
activity_ROS_M1_step      = (activity_ROS_max-activity_ROS_M1_baseline)/cc_phagocyte

epithelium  = data.frame(x = seq(1,grid_size,1),
                         y = rep(0, grid_size),
                         level_injury = 0, id = seq(1,grid_size))

epithelium[injury_site,]$level_injury = 1 # start with 1

#phagocytes
phagocyte_x = sample(1:grid_size, n_phagocytes, TRUE)
phagocyte_y = sample(2:grid_size, n_phagocytes, TRUE) 
phagocyte_pathogens_engulfed  = rep(0, n_phagocytes)
phagocyte_commensals_engulfed = rep(0, n_phagocytes)
phagocyte_num_times_activated = rep(0, n_phagocytes)
phagocyte_phenotype       = rep(0, n_phagocytes)  # 0=M0, 1=M1, 2=M2
phagocyte_activity_ROS    = rep(activity_ROS_M0_baseline, n_phagocytes)
phagocyte_activity_engulf = rep(activity_engulf_M0_baseline, n_phagocytes)
phagocyte_active_age      = rep(0, n_phagocytes)

# For bacteria_registry, use a matrix instead of list of vectors
phagocyte_bacteria_registry = matrix(0, nrow = n_phagocytes, ncol = cc_phagocyte)

# tregs
treg_x          = sample(1:grid_size, n_tregs, TRUE)
treg_y          = sample(2:grid_size, n_tregs, TRUE)
treg_active_age = rep(0, n_tregs)
treg_phenotype  = rep(0, n_tregs)  # 0=resting, 1=activated
treg_activity_SAMPs_binary = rep(0, n_tregs)

if(n_pathogens_lp == 0) {
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

# Commensals - use matrix
commensal_coords = matrix(c(
  sample(1:grid_size, n_commensals_lp, TRUE),
  sample(1:grid_size, n_commensals_lp, TRUE),
  seq(1, n_commensals_lp)
), ncol = 3)
colnames(commensal_coords) = c("x", "y", "id")

last_id_pathogen  = n_pathogens_lp
last_id_commensal = n_commensals_lp

p_prev_mic = NULL
p_prev_lym = NULL

pathogens_killed_by_ROS = 0
pathogens_killed_by_Mac = rep(0,3)

commensals_killed_by_ROS = 0
commensals_killed_by_Mac = rep(0,3)

### KEEP THE NUMBERS
epithelium_longitudinal  = matrix(0,nrow=t_max, ncol=(max_level_injury+1))
macrophages_longitudinal = matrix(0,nrow=t_max, ncol=1+2*(cc_phagocyte+1)) # M0, M1, M2 levels
microbes_longitudinal    = matrix(0,nrow=t_max, ncol=2) # first column commensals, second column pathogens
tregs_longitudinal       = matrix(0,nrow=t_max, ncol=2) # two phenotypes, active and resting
microbes_cumdeath_longitudinal = matrix(0,nrow=t_max, ncol=2*4) # first 4 columns commensals, second 4 columns pathogens - death by ROS, M0, M1, M2

for (t in 1:t_max) {
  
  # Update injury site
  injury_site_updated = which(epithelium$level_injury>0)
  
  # Update SAMPs based on activated tregs
  active_tregs = which(treg_phenotype == 1)
  if(length(active_tregs) > 0) {
    for(i in active_tregs) {
      SAMPs[treg_y[i], treg_x[i]] = SAMPs[treg_y[i], treg_x[i]] + treg_activity_SAMPs_binary[i] * add_SAMPs
    }
  }
  
  # Update ROS based on M1 phagocytes  
  M1_phagocytes = which(phagocyte_phenotype == 1)
  if(length(M1_phagocytes) > 0) {
    for(i in M1_phagocytes) {
      ROS[phagocyte_y[i], phagocyte_x[i]] = ROS[phagocyte_y[i], phagocyte_x[i]] + phagocyte_activity_ROS[i] * add_ROS
    }
  }
  
  # Move pathogens and commensals randomly (optimized for matrices)
  if(nrow(pathogen_coords) > 0) {
    # Initialize dy vector
    dy = ifelse(pathogen_coords[, "y"] == 1,
                sample(c(1), size = nrow(pathogen_coords), replace = TRUE),
                sample(c(-1, 0, 1), size = nrow(pathogen_coords), replace = TRUE))
    
    # dx sampled conditionally to avoid (0,0)
    dx = iszero_coordinates(dy)
    
    # Update positions
    pathogen_coords[, "x"] = pmin(pmax(pathogen_coords[, "x"] + dx, 1), grid_size)
    pathogen_coords[, "y"] = pmin(pmax(pathogen_coords[, "y"] + dy, 1), grid_size)
  }
  if(nrow(commensal_coords) > 0) {
    # Initialize dy vector
    dy = ifelse(commensal_coords[, "y"] == 1,
                sample(c(1), size = nrow(commensal_coords), replace = TRUE),
                sample(c(-1, 0, 1), size = nrow(commensal_coords), replace = TRUE))
    
    # dx sampled conditionally to avoid (0,0)
    dx = iszero_coordinates(dy)
    
    # Update positions
    commensal_coords[, "x"] = pmin(pmax(commensal_coords[, "x"] + dx, 1), grid_size)
    commensal_coords[, "y"] = pmin(pmax(commensal_coords[, "y"] + dy, 1), grid_size)   
  }
  
  # Pre-calculate pathogen counts touching epithelium (y=1) for each x position
  pathogen_epithelium_counts = rep(0, grid_size)
  if(nrow(pathogen_coords) > 0) {
    epithelium_pathogens = pathogen_coords[pathogen_coords[, "y"] == 1, , drop = FALSE]
    if(nrow(epithelium_pathogens) > 0) {
      pathogen_epithelium_counts = tabulate(epithelium_pathogens[, "x"], nbins = grid_size)
    }
  }
  
  # Pre-calculate commensal counts touching epithelium (y=1) for each x position
  commensal_epithelium_counts = rep(0, grid_size)
  if(nrow(commensal_coords) > 0) {
    epithelium_commensals = commensal_coords[commensal_coords[, "y"] == 1, , drop = FALSE]
    if(nrow(epithelium_commensals) > 0) {
      commensal_epithelium_counts = tabulate(epithelium_commensals[, "x"], nbins = grid_size)
    }
  }
  
  # Update DAMPs
  DAMPs[1,] = DAMPs[1,] + epithelium$level_injury*add_DAMPs
  DAMPs[1,] = DAMPs[1,] + 1*logistic_scaled_0_to_5_quantized(pathogen_epithelium_counts+commensal_epithelium_counts)*add_DAMPs
  
  # DAMPs also from pathogens? not sure. 
  DAMPs_add  = matrix(0, nrow = nrow(DAMPs), ncol = ncol(DAMPs))
  pat_counts = table(pathogen_coords[,"x"], pathogen_coords[,"y"])
  if(dim(pat_counts)[1]>0){
    pat_df = as.data.frame(pat_counts)
    names(pat_df) = c("x", "y", "count")
    pat_df$x = as.numeric(as.character(pat_df$x))
    pat_df$y = as.numeric(as.character(pat_df$y))
    pat_df$val = 0*add_DAMPs*logistic_scaled_0_to_5_quantized(pat_df$count)
    DAMPs_add[cbind(pat_df$x, pat_df$y)] = pat_df$val
  }
  DAMPs = DAMPs + DAMPs_add
  
  # Diffuse & decay DAMPs, SAMPs, and ROS.
  DAMPs   = diffuse_matrix(DAMPs, diffusion_speed_DAMPs, max_cell_value_DAMPs)
  SAMPs   = diffuse_matrix(SAMPs, diffusion_speed_SAMPs, max_cell_value_SAMPs)
  ROS     = diffuse_matrix(ROS, diffusion_speed_ROS, max_cell_value_ROS)
  DAMPs   = DAMPs - DAMPs_decay*DAMPs
  SAMPs   = SAMPs - SAMPs_decay*SAMPs
  ROS     = ROS - ros_decay*ROS
  
  if(randomize_tregs==0){
    density_matrix_tregs = DAMPs # here you can choose another density matrix if you like
  }else{
    density_matrix_tregs = 0*DAMPs # all zeros
  }
  density_matrix_phagocytes = DAMPs
  
  # Move phagocytes and tregs based on DAMPs gradient
  all_equal_tregs      = all(density_matrix_tregs == density_matrix_tregs[1, 1])
  all_equal_phagocytes = all(density_matrix_phagocytes == density_matrix_phagocytes[1, 1])
  
  ##### Tregs
  if(!all_equal_tregs){ 
    # Move tregs
    for(i in 1:length(treg_x)) {
      x = treg_x[i]
      y = treg_y[i]
      
      # Get 3x3 neighborhood (with boundary check) - same as original
      x_range = max(1, x - 1):min(grid_size, x + 1)
      y_range = max(1, y - 1):min(grid_size, y + 1)
      
      # Create all neighbor combinations
      neighbors_x = rep(x_range, each = length(y_range))
      neighbors_y = rep(y_range, times = length(x_range))
      
      # Get density values for those cells (vectorized)
      neighbor_densities = density_matrix_tregs[cbind(neighbors_y, neighbors_x)]
      
      # Normalize to get probabilities
      total = sum(neighbor_densities)
      if(total > 0) {
        probs = neighbor_densities / total
      } else {
        probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
      }
      
      # Sample a move weighted by local density
      chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
      
      # Update position (modifies in place using)
      treg_x[i] = neighbors_x[chosen_idx]
      treg_y[i] = neighbors_y[chosen_idx]
    }
  } else {
    # Random movement when no gradient 
    # tregs
    dy_treg = ifelse(treg_y == 1,
                     sample(c(1), size = length(treg_y), replace = TRUE),
                     sample(c(-1, 0, 1), size = length(treg_y), replace = TRUE))
    dx_treg = iszero_coordinates(dy_treg)# dx sampled conditionally to avoid (0,0)
    # Update positions with boundary constraints 
    treg_x = pmin(pmax(treg_x + dx_treg, 1), grid_size)
    treg_y = pmin(pmax(treg_y + dy_treg, 1), grid_size)
  }
  
  ##### Phagocytes
  if(!all_equal_phagocytes){ 
    # Move all phagocytes
    for(i in 1:length(phagocyte_x)){
      x = phagocyte_x[i]
      y = phagocyte_y[i]
      
      # Get 3x3 neighborhood (with boundary check) - same as original
      x_range = max(1, x - 1):min(grid_size, x + 1)
      y_range = max(1, y - 1):min(grid_size, y + 1)
      
      # Create all neighbor combinations
      neighbors_x = rep(x_range, each = length(y_range))
      neighbors_y = rep(y_range, times = length(x_range))
      
      # Get density values for those cells (vectorized)
      neighbor_densities = density_matrix_phagocytes[cbind(neighbors_y, neighbors_x)]
      
      # Normalize to get probabilities
      total = sum(neighbor_densities)
      if(total > 0) {
        probs = neighbor_densities / total
      } else {
        probs = rep(1 / length(neighbor_densities), length(neighbor_densities))
      }
      
      # Sample a move weighted by local density
      chosen_idx = sample(1:length(neighbors_x), 1, prob = probs)
      
      # Update position
      phagocyte_x[i] = neighbors_x[chosen_idx]
      phagocyte_y[i] = neighbors_y[chosen_idx]
    }
    
  } else {
    # Random movement when no gradient 
    # All phagocytes
    dy_phagocyte = ifelse(phagocyte_y == 1,
                          sample(c(1), size = length(phagocyte_y), replace = TRUE),
                          sample(c(-1, 0, 1), size = length(phagocyte_y), replace = TRUE))
    dx_phagocyte = iszero_coordinates(dy_phagocyte)# dx sampled conditionally to avoid (0,0)
    # Update positions with boundary constraints 
    phagocyte_x = pmin(pmax(phagocyte_x + dx_phagocyte, 1), grid_size)
    phagocyte_y = pmin(pmax(phagocyte_y + dy_phagocyte, 1), grid_size)
  }
  
  # Add new microbes based on the injured epithelium
  n_pathogens_lp_new = round(mean(epithelium$level_injury) * rate_leak_pathogen_injury * length(injury_site_updated))
  if(n_pathogens_lp_new > 0) {
    new_pathogen_coords = matrix(c(
      sample(1:grid_size, n_pathogens_lp_new, replace = TRUE, prob = epithelium$level_injury),
      rep(1, n_pathogens_lp_new),
      last_id_pathogen + seq(1, n_pathogens_lp_new)
    ), ncol = 3)
    colnames(new_pathogen_coords) = c("x", "y", "id")
    
    pathogen_coords = rbind(pathogen_coords, new_pathogen_coords)
    last_id_pathogen = last_id_pathogen + n_pathogens_lp_new
  }
  
  n_commensals_lp_new_injury   = round(mean(epithelium$level_injury) * rate_leak_commensal_injury * length(injury_site_updated))
  n_commensals_lp_new_baseline = round(rate_leak_commensal_baseline * grid_size)
  
  total_new_commensals = n_commensals_lp_new_baseline + n_commensals_lp_new_injury
  if(total_new_commensals > 0) {
    # Baseline commensals
    baseline_x = sample(1:grid_size, n_commensals_lp_new_baseline, TRUE)
    # Injury-site commensals  
    injury_x = if(n_commensals_lp_new_injury > 0) {
      sample(1:grid_size, n_commensals_lp_new_injury, TRUE, prob = epithelium$level_injury)
    } else { numeric(0) }
    
    new_commensal_coords = matrix(c(
      c(baseline_x, injury_x),
      rep(1, total_new_commensals),
      last_id_commensal + seq(1, total_new_commensals)
    ), ncol = 3)
    colnames(new_commensal_coords) = c("x", "y", "id")
    
    commensal_coords = rbind(commensal_coords, new_commensal_coords)
    last_id_commensal = last_id_commensal + total_new_commensals
  }
  
  # Update the phagocyte phenotypes
  # Vectorized operations where possible
  M0_indices = which(phagocyte_phenotype == 0)
  M1_indices = which(phagocyte_phenotype == 1)
  M2_indices = which(phagocyte_phenotype == 2)
  
  # Registry shifting every digestion_time steps (vectorized)
  if(t %% digestion_time == 0) {
    # Shift all registries at once using matrix operations
    phagocyte_bacteria_registry = cbind(
      matrix(0, nrow = nrow(phagocyte_bacteria_registry), ncol = 1),
      phagocyte_bacteria_registry[, -ncol(phagocyte_bacteria_registry)]
    )
  }
  
  # Process M0 phagocytes (candidates for activation)
  if(length(M0_indices) > 0) {
    for(i in M0_indices) {
      # Get signals (this function call is unavoidable but now faster)
      avg_DAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_DAMPs, DAMPs)
      avg_SAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_SAMPs, SAMPs)
      bacteria_count = sum(phagocyte_bacteria_registry[i, ])
      
      # Activation logic
      if(avg_DAMPs >= activation_threshold_DAMPs && avg_DAMPs > avg_SAMPs) {
        phagocyte_phenotype[i] = 1  # M1
        phagocyte_active_age[i] = 1
        phagocyte_activity_ROS[i] = activity_ROS_M1_baseline + activity_ROS_M1_step * bacteria_count
        phagocyte_activity_engulf[i] = activity_engulf_M1_baseline + activity_engulf_M1_step * bacteria_count
      } else if(avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > avg_DAMPs) {
        phagocyte_phenotype[i] = 2  # M2
        phagocyte_active_age[i] = 1
        phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
        phagocyte_activity_engulf[i] = activity_engulf_M2_baseline + activity_engulf_M2_step * bacteria_count
      }
    }
  }
  
  # Process M1/M2 phagocytes
  active_indices = c(M1_indices, M2_indices)
  if(length(active_indices) > 0) {
    # Vectorized age increment
    phagocyte_active_age[active_indices] = phagocyte_active_age[active_indices] + 1
    
    # Check which ones are old enough to potentially change
    old_enough = phagocyte_active_age[active_indices] >= active_age_limit
    candidates = active_indices[old_enough]
    
    for(i in candidates) {
      avg_DAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_DAMPs, DAMPs)
      avg_SAMPs = get_8n_avg_signal_fast(phagocyte_x[i], phagocyte_y[i], act_radius_SAMPs, SAMPs)
      
      bacteria_count = sum(phagocyte_bacteria_registry[i, ])
      
      if(avg_DAMPs >= activation_threshold_DAMPs && avg_DAMPs > avg_SAMPs) {
        phagocyte_phenotype[i] = 1
        phagocyte_active_age[i] = 1
        phagocyte_activity_ROS[i] = activity_ROS_M1_baseline + activity_ROS_M1_step * bacteria_count
        phagocyte_activity_engulf[i] = activity_engulf_M1_baseline + activity_engulf_M1_step * bacteria_count
      } else if(avg_SAMPs >= activation_threshold_SAMPs && avg_SAMPs > avg_DAMPs) {
        phagocyte_phenotype[i] = 2
        phagocyte_active_age[i] = 1
        phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
        phagocyte_activity_engulf[i] = activity_engulf_M2_baseline + activity_engulf_M2_step * bacteria_count
      } else if(avg_SAMPs < activation_threshold_SAMPs && avg_DAMPs < activation_threshold_DAMPs) {
        phagocyte_phenotype[i] = 0  # Back to M0
        phagocyte_active_age[i] = 0
        phagocyte_activity_ROS[i] = activity_ROS_M0_baseline
        phagocyte_activity_engulf[i] = activity_engulf_M0_baseline
      }
    }
  }
  
  # Treg active age updated
  active_treg_indices = which(treg_phenotype == 1)
  if(length(active_treg_indices) > 0) {
    # Vectorized age check
    old_tregs = active_treg_indices[treg_active_age[active_treg_indices] >= active_age_limit]
    young_tregs = active_treg_indices[treg_active_age[active_treg_indices] < active_age_limit]
    
    # Age the young ones
    if(length(young_tregs) > 0) {
      treg_active_age[young_tregs] = treg_active_age[young_tregs] + 1
    }
    
    # Shut down the old ones
    if(length(old_tregs) > 0) {
      treg_phenotype[old_tregs] = 0
      treg_active_age[old_tregs] = 0
      treg_activity_SAMPs_binary[old_tregs] = 0
    }
  }
  
  # Engulfment process
  # Pre-calculate all phagocyte positions for faster lookup
  phagocyte_positions = paste(phagocyte_x, phagocyte_y, sep = "_")
  
  for(i in 1:length(phagocyte_x)) {
    px = phagocyte_x[i]
    py = phagocyte_y[i]
    
    # Fast pathogen overlap check using vectorized operations
    if(nrow(pathogen_coords) > 0) {
      pathogen_overlap = (pathogen_coords[, "x"] == px) & (pathogen_coords[, "y"] == py)
      pathogen_indices = which(pathogen_overlap)
      
      if(length(pathogen_indices) > 0) {
        # Stochastic engulfment
        engulf_success = runif(length(pathogen_indices)) < phagocyte_activity_engulf[i]
        indices_to_engulf = pathogen_indices[engulf_success]
        
        if(length(indices_to_engulf) > 0) {
          phagocyte_pathogens_engulfed[i] = phagocyte_pathogens_engulfed[i] + length(indices_to_engulf)
          
          # Remove engulfed pathogens (create new matrix without these rows)
          pathogen_coords = pathogen_coords[-indices_to_engulf, , drop = FALSE]
          
          # Update bacteria registry efficiently
          phagocyte_bacteria_registry[i, ] = shift_insert_fast(phagocyte_bacteria_registry[i, ], 
                                                               rep(1, length(indices_to_engulf)))
          
          # Update kill count
          phagocyte_phenotype_index = phagocyte_phenotype[i] + 1
          pathogens_killed_by_Mac[phagocyte_phenotype_index] = pathogens_killed_by_Mac[phagocyte_phenotype_index] + length(indices_to_engulf)
        }
      }
    }
    
    # Similar logic for commensals
    if(nrow(commensal_coords) > 0) {
      commensal_overlap = (commensal_coords[, "x"] == px) & (commensal_coords[, "y"] == py)
      commensal_indices = which(commensal_overlap)
      
      if(length(commensal_indices) > 0) {
        engulf_success = runif(length(commensal_indices)) < phagocyte_activity_engulf[i]
        indices_to_engulf = commensal_indices[engulf_success]
        
        if(length(indices_to_engulf) > 0) {
          phagocyte_commensals_engulfed[i] = phagocyte_commensals_engulfed[i] + length(indices_to_engulf)
          commensal_coords = commensal_coords[-indices_to_engulf, , drop = FALSE]
          
          phagocyte_bacteria_registry[i, ] = shift_insert_fast(phagocyte_bacteria_registry[i, ], 
                                                               rep(1, length(indices_to_engulf)))
          
          phagocyte_phenotype_index = phagocyte_phenotype[i] + 1
          commensals_killed_by_Mac[phagocyte_phenotype_index] = commensals_killed_by_Mac[phagocyte_phenotype_index] + length(indices_to_engulf)
        }
      }
    }
  }
  
  # Treg activation & effector actions
  if(allow_tregs_to_do_their_job) {
    # Get all M1 phagocytes for treg interactions
    M1_phagocyte_indices = which(phagocyte_phenotype == 1)
    
    if(length(M1_phagocyte_indices) > 0) {
      for(i in M1_phagocyte_indices) {
        px = phagocyte_x[i]
        py = phagocyte_y[i]
        
        # Check if a Treg is in the vicinity (vectorized)
        treg_distances_x = abs(treg_x - px)
        treg_distances_y = abs(treg_y - py)
        nearby_treg_indices = which(treg_distances_x <= treg_vicinity_effect & 
                                      treg_distances_y <= treg_vicinity_effect)
        
        if(length(nearby_treg_indices) > 0) {
          num_pat_antigens = phagocyte_pathogens_engulfed[i]
          num_com_antigens = phagocyte_commensals_engulfed[i]
          
          if((num_pat_antigens + num_com_antigens) > 0) {
            
            # com_pat_vec = rep(0,num_pat_antigens + num_com_antigens)
            # com_pat_vec[1:num_com_antigens]=1 # commensals are 1, pathogens are 0
            # pat_indices = which(com_pat_vec == 0)
            # com_pat_vec_flipped = com_pat_vec
            # if (length(pat_indices) > 0) {
            #   flip_pathogens = runif(length(pat_indices)) < (1 - treg_discrimination_efficiency)
            #   com_pat_vec_flipped[pat_indices[flip_pathogens]] = 1
            # }
            # 
            # # rat_com_pat = num_com_antigens / (num_pat_antigens + num_com_antigens)
            # rat_com_pat = mean(com_pat_vec_flipped)
            
            rat_com_pat_real = num_com_antigens/(num_com_antigens+num_pat_antigens)
            alpha = (1-treg_discrimination_efficiency)*1 + treg_discrimination_efficiency*(rat_com_pat_real*precision)
            beta  = (1-treg_discrimination_efficiency)*1 + treg_discrimination_efficiency*((1 - rat_com_pat_real)*precision)
            # rat_com_pat = rbeta(1, alpha, beta) # built in
            rat_com_pat = sample_rbeta(alpha, beta) # to match python
            
            if(rat_com_pat > rat_com_pat_threshold) {
              # TREG ACTIVATION (vectorized for all nearby tregs)
              treg_phenotype[nearby_treg_indices] = 1
              treg_activity_SAMPs_binary[nearby_treg_indices] = 1
              treg_active_age[nearby_treg_indices] = 1
              
              if(allow_tregs_to_suppress_cognate) {
                # Suppress this M1 phagocyte
                phagocyte_phenotype[i] = 2  # Convert to M2
                phagocyte_active_age[i] = 1
                bacteria_count = sum(phagocyte_bacteria_registry[i, ])
                phagocyte_activity_ROS[i] = activity_ROS_M2_baseline
                phagocyte_activity_engulf[i] = activity_engulf_M2_baseline + 
                  activity_engulf_M2_step * bacteria_count
              }
            }
          }
        }
      }
    }
  }
  
  # Kill pathogens with ROS
  if(nrow(pathogen_coords) > 0) {
    # Vectorized ROS exposure calculation
    pathogen_avg_ROS = numeric(nrow(pathogen_coords))
    
    for(i in 1:nrow(pathogen_coords)) {
      pathogen_avg_ROS[i] = get_8n_avg_signal_fast(pathogen_coords[i, "x"], 
                                                   pathogen_coords[i, "y"], 
                                                   act_radius_ROS, ROS)
    }
    
    # Vectorized death determination
    pathogens_to_kill = which(pathogen_avg_ROS > th_ROS_microbe)
    
    if(length(pathogens_to_kill) > 0) {
      pathogen_coords = pathogen_coords[-pathogens_to_kill, , drop = FALSE]
      pathogens_killed_by_ROS = pathogens_killed_by_ROS + length(pathogens_to_kill)
    }
  }
  
  # Kill commensals with ROS
  if(nrow(commensal_coords) > 0) {
    # Vectorized ROS exposure calculation
    commensal_avg_ROS = numeric(nrow(commensal_coords))
    
    for(i in 1:nrow(commensal_coords)) {
      commensal_avg_ROS[i] = get_8n_avg_signal_fast(commensal_coords[i, "x"], 
                                                    commensal_coords[i, "y"], 
                                                    act_radius_ROS, ROS)
    }
    
    # Vectorized death determination
    commensals_to_kill = which(commensal_avg_ROS > th_ROS_microbe)
    
    if(length(commensals_to_kill) > 0) {
      commensal_coords = commensal_coords[-commensals_to_kill, , drop = FALSE]
      commensals_killed_by_ROS = commensals_killed_by_ROS + length(commensals_to_kill)
    }
  }
  
  # Injur epithelium with pathogens / commensals
  for(i in 1:nrow(epithelium)) {
    px = epithelium$x[i]
    
    # Get ROS values in vicinity
    x_coordinates = pmax(1, pmin(grid_size, (px - act_radius_ROS):(px + act_radius_ROS)))
    ros_values = ROS[1, x_coordinates]  # Row right below epithelium
    mean_ros = mean(ros_values)
    
    # Increase level_injury based on pathogen count
    count_pathogens = pathogen_epithelium_counts[px]
    # epithelium$level_injury[i] = epithelium$level_injury[i] + round(log(count_pathogens + 1))
    epithelium$level_injury[i] = epithelium$level_injury[i] + logistic_scaled_0_to_5_quantized(count_pathogens)
    
    # RULE 2: Increase level_injury based on ROS
    if(mean_ros > th_ROS_epith_recover) {
      epithelium$level_injury[i] = epithelium$level_injury[i] + 1
    }
    
    # Apply maximum injury constraint
    epithelium$level_injury[i] = min(epithelium$level_injury[i], max_level_injury)
    
    # RECOVERY: Stochastic recovery when injured
    if(epithelium$level_injury[i] > 0 && runif(1) < epith_recovery_chance) {
      epithelium$level_injury[i] = max(0, epithelium$level_injury[i] - 1)
    }
  }
  
  # Save abundances
  epithelium_longitudinal[t, ] = as.numeric(table(factor(epithelium$level_injury, levels = 0:5)))
  
  # Phenotype counting
  phagocyte_counts = c(
    sum(phagocyte_phenotype == 0),  # M0
    tabulate(phagocyte_active_age[phagocyte_phenotype == 1] + 1, cc_phagocyte+1),  # M1 by level
    tabulate(phagocyte_active_age[phagocyte_phenotype == 2] + 1, cc_phagocyte+1)   # M2 by level  
  )
  macrophages_longitudinal[t, ] = phagocyte_counts
  
  microbes_longitudinal[t, ] = c(nrow(commensal_coords), nrow(pathogen_coords))
  tregs_longitudinal[t, ] = c(sum(treg_phenotype == 0), sum(treg_phenotype == 1))
  microbes_cumdeath_longitudinal[t, ] = c(commensals_killed_by_ROS, commensals_killed_by_Mac, 
                                          pathogens_killed_by_ROS, pathogens_killed_by_Mac)
}

# Combine all matrices into one
longitudinal_df = data.frame(
  epithelium_longitudinal,
  macrophages_longitudinal,
  microbes_longitudinal,
  tregs_longitudinal,
  microbes_cumdeath_longitudinal
)

# Add time column
longitudinal_df$t = 1:t_max
# # Add metadata columns (same value repeated for each row)
longitudinal_df$sterile                       = sterile
longitudinal_df$allow_tregs_to_do_their_job   = allow_tregs_to_do_their_job
longitudinal_df$allow_tregs_to_suppress_cognate = allow_tregs_to_suppress_cognate
longitudinal_df$randomize_tregs                 = randomize_tregs
#--------------------------------------------------------------------------------------------------------
# Optional: Reorder columns to have metadata first
longitudinal_df = longitudinal_df %>% dplyr::select(t,sterile,allow_tregs_to_do_their_job,
                                                    allow_tregs_to_suppress_cognate,
                                                    randomize_tregs,
                                                    everything())

colnames_insert = c('epithelial_healthy','epithelial_inj_1','epithelial_inj_2','epithelial_inj_3','epithelial_inj_4','epithelial_inj_5',
                    'phagocyte_M0','phagocyte_M1_L_0','phagocyte_M1_L_1','phagocyte_M1_L_2','phagocyte_M1_L_3','phagocyte_M1_L_4','phagocyte_M1_L_5',
                    'phagocyte_M2_L_0','phagocyte_M2_L_1','phagocyte_M2_L_2','phagocyte_M2_L_3','phagocyte_M2_L_4','phagocyte_M2_L_5',
                    'commensal','pathogen','treg_resting','treg_active','C_ROS','C_M0','C_M1','C_M2','P_ROS','P_M0','P_M1','P_M2')

colnames(longitudinal_df)[6:36] = colnames_insert #+30

# ------ PLOTTING FUNCTIONS ------------------------------------------------------------------------  
plot_faceted_8 = function(data, variables, title) {
  data_long = data %>%
    dplyr::select(t, sterile, allow_tregs_to_do_their_job, randomize_tregs, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")
  
  p = ggplot(data_long, aes(x = t, y = value, color = variable)) +
    geom_line(alpha = 1, linewidth = 1) +
    # facet_grid(sterile ~ allow_tregs_to_do_their_job + randomize_tregs, labeller = label_both) +
    facet_grid(randomize_tregs ~ allow_tregs_to_do_their_job + sterile, labeller = label_both) +
    scale_color_manual(values = agent_colors) +
    theme_minimal() +
    labs(title = title, x = "Time", y = "Count", color = "Agent")
  
  return(p)
}

agent_colors = c(
  epithelial_healthy = "#B0E2FF",
  epithelial_inj_1   = "#8CB4E5",
  epithelial_inj_2   = "#6987CC",
  epithelial_inj_3   = "#465AB2",
  epithelial_inj_4   = "#232D99",
  epithelial_inj_5   = "#000080",
  phagocyte_M0       = "grey70",
  phagocyte_M1_L_0   = "#F8C8E8",
  phagocyte_M1_L_1   = "#F397D6",
  phagocyte_M1_L_2   = "#E754C4",
  phagocyte_M1_L_3   = "#D12CA0",
  phagocyte_M1_L_4   = "#A5177A",
  phagocyte_M1_L_5   = "#6B0C4F",
  phagocyte_M2_L_0   = "#CDEFE3", 
  phagocyte_M2_L_1   = "#97D6BC",  
  phagocyte_M2_L_2   = "#61BD96",  
  phagocyte_M2_L_3   = "#3BA578",  
  phagocyte_M2_L_4   = "#2E8B57",  
  phagocyte_M2_L_5   = "#1F5C3B",  
  phagocyte_M0   = "grey70",
  phagocyte_M1   = "#E754C4",
  phagocyte_M2   = "#3BA578",
  treg_resting = "#D8BFD8",
  treg_active  = "#967BB6",
  commensal    = "turquoise2",
  pathogen     = "firebrick1",
  C_ROS = "black",
  C_M0  = "grey70",
  C_M1  = "#D12CA0",
  C_M2  = "#3BA578",
  P_ROS = "black",
  P_M0  = "grey70",
  P_M1  = "#D12CA0",
  P_M2  = "#3BA578"
)
# ------ PLOTTING FUNCTIONS ------------------------------------------------------------------------  

# ------ PLOTTING ------------------------------------------------------------------------  

p_epithelium = plot_faceted_8(longitudinal_df, c("epithelial_healthy", paste0("epithelial_inj_", 1:5)),"Epithelial Cell Dynamics")
print(p_epithelium)

p_microbes= plot_faceted_8(longitudinal_df, c("commensal", "pathogen"), "Microbe Dynamics")
print(p_microbes)

p_tregs= plot_faceted_8(longitudinal_df, c("treg_resting", "treg_active"), "Treg Dynamics")
print(p_tregs)

longitudinal_df = longitudinal_df %>% dplyr::mutate(phagocyte_M1 = phagocyte_M1_L_0+phagocyte_M1_L_1+phagocyte_M1_L_2+phagocyte_M1_L_3+phagocyte_M1_L_4+phagocyte_M1_L_5)
longitudinal_df = longitudinal_df %>% dplyr::mutate(phagocyte_M2 = phagocyte_M2_L_0+phagocyte_M2_L_1+phagocyte_M2_L_2+phagocyte_M2_L_3+phagocyte_M2_L_4+phagocyte_M2_L_5)
p_phagocytes= plot_faceted_8(longitudinal_df, c("phagocyte_M0", "phagocyte_M1", "phagocyte_M2"), "Phagocyte Dynamics")
print(p_phagocytes)


###

p_results = read.csv('/Users/burcutepekule/Dropbox/Treg_problem_v2/python/simulation_results.csv')
p_results = p_results %>% dplyr::mutate(t=t+1)

library(cowplot)

# Basic 2x2 grid
p_all=plot_grid(p_epithelium, 
          p_microbes, 
          p_phagocytes, 
          p_tregs, 
          ncol = 2, 
          nrow = 2)
ggsave("/Users/burcutepekule/Desktop/tregs/simulation_results_R.png", 
       plot = p_all,
       width = 10,     # width in inches
       height = 8,     # height in inches
       dpi = 300)      # resolution (dots per inch)
