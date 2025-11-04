rm(list=ls())

# Load required libraries
library(keras)
library(tensorflow)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(purrr)

# ===============================================
# 1. DATA PREPROCESSING FUNCTIONS
# ===============================================

prepare_epithelial_score <- function(df) {
  # Calculate epithelial score for each time point
  df %>%
    mutate(epithelial_score = 6 * epithelial_healthy + 
             5 * epithelial_inj_1 + 
             4 * epithelial_inj_2 + 
             3 * epithelial_inj_3 + 
             2 * epithelial_inj_4 + 
             1 * epithelial_inj_5)
}

create_scenario_encoding <- function(sterile, allow_tregs, randomize_tregs) {
  # Create one-hot encoding for scenarios
  scenario_id <- sterile * 4 + allow_tregs * 2 + randomize_tregs
  scenario_onehot <- matrix(0, nrow = length(scenario_id), ncol = 6)
  for(i in 1:length(scenario_id)) {
    if(allow_tregs[i] == 0) {
      # Only 2 scenarios when tregs are off (sterile or not)
      scenario_onehot[i, sterile[i] + 1] <- 1
    } else {
      # 4 scenarios when tregs are on
      scenario_onehot[i, scenario_id[i] + 1] <- 1
    }
  }
  scenario_onehot
}

load_and_process_batch <- function(param_set_ids, base_path, n_time_points = 500) {
  
  batch_data <- list()
  
  for(param_id in param_set_ids) {
    file_path <- paste0(base_path, '/simulation_results_param_set_', param_id, '.csv')
    
    if(!file.exists(file_path)) next
    
    # Load results
    results <- read_csv(file_path, show_col_types = FALSE) %>%
      left_join(parameters_df, by = 'param_set_id')
    
    # Calculate epithelial scores
    results <- prepare_epithelial_score(results)
    
    # Get unique scenarios for this parameter set
    scenarios <- results %>%
      select(sterile, allow_tregs_to_do_their_job, randomize_tregs) %>%
      rename(allow_tregs = allow_tregs_to_do_their_job) %>%
      distinct()
    
    for(i in 1:nrow(scenarios)) {
      scenario_data <- results %>%
        filter(sterile == scenarios$sterile[i],
               allow_tregs_to_do_their_job == scenarios$allow_tregs[i],
               randomize_tregs == scenarios$randomize_tregs[i])
      
      # Extract time series (assuming 500 time points)
      time_series <- scenario_data$epithelial_score[1:n_time_points]
      
      # Extract parameters (first row since they're constant)
      params <- scenario_data[1, parameter_names] %>% as.numeric()
      
      # Create scenario encoding
      scenario_enc <- create_scenario_encoding(
        scenarios$sterile[i], 
        scenarios$allow_tregs[i], 
        scenarios$randomize_tregs[i]
      )
      
      batch_data[[length(batch_data) + 1]] <- list(
        param_set_id = param_id,
        time_series = time_series,
        parameters = params,
        scenario = scenario_enc,
        sterile = scenarios$sterile[i],
        allow_tregs = scenarios$allow_tregs[i],
        randomize_tregs = scenarios$randomize_tregs[i]
      )
    }
  }
  
  batch_data
}

# ===============================================
# 2. CONDITIONAL VAE MODEL IN R/KERAS
# ===============================================


# Working version for Keras 2.15.0
create_simple_vae <- function(time_steps = 500, 
                              n_params = 24, 
                              n_scenarios = 6, 
                              latent_dim = 32) {
  
  # Clear any previous models
  k_clear_session()
  
  # ENCODER
  time_series_input <- layer_input(shape = c(time_steps), name = "time_series")
  params_input <- layer_input(shape = c(n_params), name = "params")
  scenario_input <- layer_input(shape = c(n_scenarios), name = "scenario")
  
  # Simple dense encoder
  encoded <- layer_concatenate(list(
    time_series_input,
    params_input,
    scenario_input
  )) %>%
    layer_dense(units = 256, activation = "relu") %>%
    layer_dense(units = 128, activation = "relu")
  
  z_mean <- encoded %>% layer_dense(units = latent_dim, name = "z_mean")
  z_log_var <- encoded %>% layer_dense(units = latent_dim, name = "z_log_var")
  
  # Sampling layer using keras functional API
  z <- layer_lambda(
    f = function(x) {
      z_mean <- x[[1]]
      z_log_var <- x[[2]]
      batch <- k_shape(z_mean)[1]
      dim <- k_int_shape(z_mean)[2]
      epsilon <- k_random_normal(shape = c(batch, dim))
      z_mean + k_exp(0.5 * z_log_var) * epsilon
    }
  )(list(z_mean, z_log_var))
  
  # DECODER
  decoded <- layer_concatenate(list(z, params_input, scenario_input)) %>%
    layer_dense(units = 128, activation = "relu") %>%
    layer_dense(units = 256, activation = "relu") %>%
    layer_dense(units = time_steps, activation = "linear", name = "output")
  
  # Create encoder model
  encoder <- keras_model(
    inputs = list(time_series_input, params_input, scenario_input),
    outputs = list(z_mean, z_log_var, z),
    name = "encoder"
  )
  
  # Create VAE model
  vae <- keras_model(
    inputs = list(time_series_input, params_input, scenario_input),
    outputs = decoded,
    name = "vae"
  )
  
  # Custom VAE loss
  vae_loss <- custom_metric("vae_loss", function(y_true, y_pred) {
    # Reconstruction loss
    reconstruction_loss <- k_mean(k_square(y_true - y_pred))
    
    # KL divergence loss
    kl_loss <- -0.5 * k_mean(
      1 + z_log_var - k_square(z_mean) - k_exp(z_log_var),
      axis = -1L
    )
    
    # Total loss
    reconstruction_loss * time_steps + kl_loss
  })
  
  # Compile using the correct syntax for Keras 2.15.0
  vae %>% keras::compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    loss = 'mse',  # Use standard MSE first
    metrics = list('mae')
  )
  
  list(vae = vae, encoder = encoder)
}
# ===============================================
# 3. TRAINING DATA GENERATOR
# ===============================================

create_data_generator <- function(param_set_ids, base_path, batch_size = 32, 
                                  n_time_points = 500, parameters_df) {
  
  function() {
    # Sample batch of parameter sets
    batch_params <- sample(param_set_ids, min(batch_size, length(param_set_ids)), replace = TRUE)
    
    # Load batch data
    batch_data <- load_and_process_batch(batch_params, base_path, n_time_points)
    
    if(length(batch_data) == 0) return(NULL)
    
    # Prepare arrays
    n_samples <- length(batch_data)
    time_series_array <- matrix(0, nrow = n_samples, ncol = n_time_points)
    params_array <- matrix(0, nrow = n_samples, ncol = 24)
    scenario_array <- matrix(0, nrow = n_samples, ncol = 6)
    
    for(i in 1:n_samples) {
      time_series_array[i,] <- batch_data[[i]]$time_series
      params_array[i,] <- batch_data[[i]]$parameters
      scenario_array[i,] <- batch_data[[i]]$scenario
    }
    
    # Normalize
    time_series_array <- scale(time_series_array)
    params_array <- scale(params_array)
    
    list(
      list(time_series_array, params_array, scenario_array),
      time_series_array  # y = x for autoencoder
    )
  }
}

# ===============================================
# 4. ANALYSIS FUNCTIONS
# ===============================================

analyze_treg_effects <- function(encoder, data_batch) {
  # Extract latent representations
  latent_reps <- encoder %>% predict(
    list(data_batch$time_series, data_batch$params, data_batch$scenarios)
  )
  
  z_mean <- latent_reps[[1]]
  
  # Compare with/without Tregs for same parameters
  treg_effects <- list()
  
  unique_params <- unique(data_batch$param_set_id)
  
  for(param_id in unique_params) {
    idx <- which(data_batch$param_set_id == param_id)
    
    if(length(idx) < 2) next
    
    # Find pairs with/without tregs
    no_treg_idx <- idx[data_batch$allow_tregs[idx] == 0]
    with_treg_idx <- idx[data_batch$allow_tregs[idx] == 1]
    
    if(length(no_treg_idx) > 0 && length(with_treg_idx) > 0) {
      # Calculate distance in latent space
      for(i in no_treg_idx) {
        for(j in with_treg_idx) {
          if(data_batch$sterile[i] == data_batch$sterile[j]) {
            latent_distance <- sqrt(sum((z_mean[i,] - z_mean[j,])^2))
            
            treg_effects[[length(treg_effects) + 1]] <- data.frame(
              param_set_id = param_id,
              sterile = data_batch$sterile[i],
              randomized = data_batch$randomize_tregs[j],
              latent_distance = latent_distance,
              parameters = t(data_batch$parameters[i,])
            )
          }
        }
      }
    }
  }
  
  bind_rows(treg_effects)
}

# ===============================================
# 5. MAIN TRAINING SCRIPT
# ===============================================

# Load parameters
parameters_df <- read_csv('/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv', 
                          show_col_types = FALSE)

base_path <- '/Users/burcutepekule/Desktop/tregs/mass_sim_results'

# Get all parameter set IDs
all_param_ids <- unique(parameters_df$param_set_id)

# Split into train/validation
n_params <- length(all_param_ids)
train_ids <- sample(all_param_ids, floor(0.8 * n_params))
val_ids <- setdiff(all_param_ids, train_ids)

# Create model
models <- create_simple_vae(
  time_steps = 500,
  n_params = 24,
  n_scenarios = 6,
  latent_dim = 32
)

vae <- models$vae
encoder <- models$encoder

# Create data generators
train_gen <- create_data_generator(train_ids, base_path, batch_size = 32, 
                                   parameters_df = parameters_df)
val_gen <- create_data_generator(val_ids, base_path, batch_size = 32,
                                 parameters_df = parameters_df)

# Train model (simplified - you'd want callbacks, early stopping, etc.)
history <- vae %>% fit(
  x = train_gen,
  steps_per_epoch = 100,
  epochs = 50,
  validation_data = val_gen,
  validation_steps = 20
)

# ===============================================
# 6. IDENTIFY BENEFICIAL TREG PARAMETERS
# ===============================================

# Load a batch for analysis
analysis_batch <- load_and_process_batch(
  sample(all_param_ids, 100), 
  base_path,
  n_time_points = 500
)

# Analyze Treg effects
treg_analysis <- analyze_treg_effects(encoder, analysis_batch)

# Find parameters where Tregs are beneficial (large latent distance = different dynamics)
beneficial_params <- treg_analysis %>%
  group_by(param_set_id, sterile) %>%
  summarize(
    mean_effect = mean(latent_distance),
    randomization_impact = latent_distance[randomized == 1] - latent_distance[randomized == 0]
  ) %>%
  arrange(desc(mean_effect))

print(beneficial_params)

# Visualize latent space
# You can use UMAP or t-SNE on the latent representations to visualize clusters