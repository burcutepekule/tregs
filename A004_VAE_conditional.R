# Clean VAE Implementation for Treg Analysis using Torch
# ======================================================

rm(list=ls())
library(torch)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# Set seed for reproducibility
set.seed(42)
torch_manual_seed(42)

# ======================================================
# 1. DATA LOADING AND PREPROCESSING
# ======================================================

# Load parameter names
parameter_names = c(
  "th_ROS_microbe", "th_ROS_epith_recover", "epith_recovery_chance",
  "rat_com_pat_threshold", "diffusion_speed_DAMPs", "diffusion_speed_SAMPs",
  "diffusion_speed_ROS", "add_ROS", "add_DAMPs", "add_SAMPs",
  "ros_decay", "DAMPs_decay", "SAMPs_decay", "activation_threshold_DAMPs",
  "activation_threshold_SAMPs", "activity_engulf_M0_baseline",
  "activity_engulf_M1_baseline", "activity_engulf_M2_baseline",
  "activity_ROS_M1_baseline", "rate_leak_commensal_injury",
  "rate_leak_pathogen_injury", "rate_leak_commensal_baseline",
  "active_age_limit", "treg_discrimination_efficiency"
)

# Function to calculate epithelial score
calculate_epithelial_score = function(df) {
  df %>%
    mutate(epithelial_score = 6 * epithelial_healthy + 
             5 * epithelial_inj_1 + 
             4 * epithelial_inj_2 + 
             3 * epithelial_inj_3 + 
             2 * epithelial_inj_4 + 
             1 * epithelial_inj_5)
}

# Function to load and process data for multiple parameter sets
load_simulation_data = function(param_set_ids, base_path, parameters_df, 
                                n_time_points = 500, n_replicates = 10) {
  
  all_data = list()
  
  for(param_id in param_set_ids) {
    file_path = paste0(base_path, '/simulation_results_param_set_', param_id, '.csv')
    results   = read_csv(file_path, show_col_types = FALSE)
    # Calculate epithelial scores
    results   = calculate_epithelial_score(results)
    
    # Get parameters for this set
    params = parameters_df %>% 
      filter(param_set_id == param_id) %>%
      dplyr::select(all_of(parameter_names))
    
    # Process each scenario
    scenarios = results %>%
      dplyr::select(sterile, allow_tregs_to_do_their_job, randomize_tregs) %>%
      distinct()
    
    for(i in 1:nrow(scenarios)) {
      for(rep in 1:n_replicates) {
        rep_data = results %>%
          filter(replicate_id == (rep-1),
                 sterile == scenarios$sterile[i],
                 allow_tregs_to_do_their_job == scenarios$allow_tregs_to_do_their_job[i],
                 randomize_tregs == scenarios$randomize_tregs[i]
          )
        
        ts_data = rep_data$epithelial_score[1:n_time_points]
        
        # Create scenario encoding (one-hot)
        scenario_vec = rep(0, 6)
        # Encoding: [no_treg_nonsterile, no_treg_sterile, treg_nonsterile_nonrandom, 
        #            treg_nonsterile_random, treg_sterile_nonrandom, treg_sterile_random]
        if(scenarios$allow_tregs_to_do_their_job[i] == 0) {
          scenario_vec[scenarios$sterile[i] + 1] = 1
        } else {
          idx = 3 + scenarios$sterile[i] * 2 + scenarios$randomize_tregs[i]
          scenario_vec[idx] = 1
        }
        
        all_data[[length(all_data) + 1]] = list(
          time_series = ts_data,
          params = as.numeric(params[1,]),
          scenario = scenario_vec,
          param_id = param_id,
          replicate_id = rep,
          sterile = scenarios$sterile[i],
          allow_tregs = scenarios$allow_tregs_to_do_their_job[i],
          randomize_tregs = scenarios$randomize_tregs[i]
        )
      }
    }
  }
  all_data
}

# Convert list data to tensors - UPDATED to include replicate info
prepare_tensors = function(data_list, normalize = TRUE) {
  n_samples = length(data_list)
  if(n_samples == 0) return(NULL)
  
  n_time = length(data_list[[1]]$time_series)
  n_params = length(data_list[[1]]$params)
  n_scenarios = length(data_list[[1]]$scenario)
  
  # Initialize matrices
  time_series_mat = matrix(0, nrow = n_samples, ncol = n_time)
  params_mat = matrix(0, nrow = n_samples, ncol = n_params)
  scenario_mat = matrix(0, nrow = n_samples, ncol = n_scenarios)
  metadata = data.frame()
  
  for(i in 1:n_samples) {
    time_series_mat[i,] = data_list[[i]]$time_series
    params_mat[i,] = data_list[[i]]$params
    scenario_mat[i,] = data_list[[i]]$scenario
    
    metadata = rbind(metadata, data.frame(
      param_id = data_list[[i]]$param_id,
      replicate_id = data_list[[i]]$replicate_id,  # ADD THIS
      sterile = data_list[[i]]$sterile,
      allow_tregs = data_list[[i]]$allow_tregs,
      randomize_tregs = data_list[[i]]$randomize_tregs
    ))
  }
  
  # Normalize if requested
  if(normalize) {
    time_series_mat = scale(time_series_mat)
    params_mat = scale(params_mat)
    
    attr(time_series_mat, "scaled:center") -> ts_center
    attr(time_series_mat, "scaled:scale") -> ts_scale
    attr(params_mat, "scaled:center") -> param_center
    attr(params_mat, "scaled:scale") -> param_scale
  }
  
  # Convert to torch tensors
  list(
    time_series = torch_tensor(time_series_mat, dtype = torch_float32()),
    params = torch_tensor(params_mat, dtype = torch_float32()),
    scenarios = torch_tensor(scenario_mat, dtype = torch_float32()),
    metadata = metadata,
    normalization = if(normalize) list(
      ts_center = ts_center, ts_scale = ts_scale,
      param_center = param_center, param_scale = param_scale
    ) else NULL
  )
}

# ======================================================
# 2. VAE MODEL DEFINITION
# ======================================================

VAE = nn_module(
  "VAE",
  
  initialize = function(time_steps = 500, n_params = 24, n_scenarios = 6, 
                        latent_dim = 32, hidden_dim = 256) {
    self$time_steps = time_steps
    self$n_params = n_params
    self$n_scenarios = n_scenarios
    self$latent_dim = latent_dim
    
    input_dim = time_steps + n_params + n_scenarios
    
    # Encoder network
    self$encoder = nn_sequential(
      nn_linear(input_dim, hidden_dim),
      nn_relu(),
      nn_dropout(0.2),
      nn_linear(hidden_dim, hidden_dim/2),
      nn_relu(),
      nn_dropout(0.2)
    )
    
    # Latent space parameters
    self$fc_mu = nn_linear(hidden_dim/2, latent_dim)
    self$fc_logvar = nn_linear(hidden_dim/2, latent_dim)
    
    # Decoder network
    decoder_input_dim = latent_dim + n_params + n_scenarios
    self$decoder = nn_sequential(
      nn_linear(decoder_input_dim, hidden_dim/2),
      nn_relu(),
      nn_dropout(0.2),
      nn_linear(hidden_dim/2, hidden_dim),
      nn_relu(),
      nn_linear(hidden_dim, time_steps)
    )
  },
  
  encode = function(x) {
    h = self$encoder(x)
    list(
      mu = self$fc_mu(h),
      logvar = self$fc_logvar(h)
    )
  },
  
  reparameterize = function(mu, logvar) {
    if(self$training) {
      std = torch_exp(0.5 * logvar)
      eps = torch_randn_like(std)
      mu + eps * std
    } else {
      mu  # Use mean during evaluation
    }
  },
  
  forward = function(time_series, params, scenarios) {
    # Concatenate all inputs
    x = torch_cat(list(time_series, params, scenarios), dim = 2)
    
    # Encode to latent space
    encoding = self$encode(x)
    mu = encoding$mu
    logvar = encoding$logvar
    
    # Sample from latent space
    z = self$reparameterize(mu, logvar)
    
    # Decode with conditioning
    decoder_input = torch_cat(list(z, params, scenarios), dim = 2)
    reconstruction = self$decoder(decoder_input)
    
    list(
      reconstruction = reconstruction,
      mu = mu,
      logvar = logvar,
      z = z
    )
  }
)

# Loss function for VAE
vae_loss = function(reconstruction, target, mu, logvar, kl_weight = 0.001) {
  # Reconstruction loss (MSE)
  recon_loss = nn_mse_loss()(reconstruction, target)
  
  # KL divergence loss
  kl_loss = -0.5 * torch_mean(1 + logvar - mu$pow(2) - logvar$exp())
  
  # Total loss
  total_loss = recon_loss + kl_weight * kl_loss
  
  list(
    total = total_loss,
    reconstruction = recon_loss,
    kl = kl_loss
  )
}

# ======================================================
# 3. TRAINING FUNCTIONS
# ======================================================

train_vae_with_validation = function(model, data_tensors, 
                                      val_split = 0.2,
                                      epochs = 100, 
                                      batch_size = 32, 
                                      learning_rate = 0.001, 
                                      kl_weight = 0.001,
                                      early_stopping_patience = 10) {
  
  optimizer = optim_adam(model$parameters, lr = learning_rate)
  
  # Split data into train/validation
  n_samples = data_tensors$time_series$shape[1]
  n_val = floor(n_samples * val_split)
  n_train = n_samples - n_val
  
  # Random split
  indices = sample(n_samples)
  train_indices = indices[1:n_train]
  val_indices = indices[(n_train+1):n_samples]
  
  # Prepare validation data
  val_ts = data_tensors$time_series[val_indices, ]
  val_params = data_tensors$params[val_indices, ]
  val_scenarios = data_tensors$scenarios[val_indices, ]
  
  history = list(
    train_loss = numeric(epochs),
    val_loss = numeric(epochs),
    train_recon = numeric(epochs),
    val_recon = numeric(epochs),
    train_kl = numeric(epochs),
    val_kl = numeric(epochs)
  )
  
  best_val_loss = Inf
  patience_counter = 0
  
  for(epoch in 1:epochs) {
    # Training phase
    model$train()
    train_losses = list(total = 0, recon = 0, kl = 0)
    
    # Shuffle training indices
    train_indices_shuffled = sample(train_indices)
    n_batches = ceiling(n_train / batch_size)
    
    for(batch in 1:n_batches) {
      start_idx = (batch - 1) * batch_size + 1
      end_idx = min(batch * batch_size, n_train)
      batch_idx = train_indices_shuffled[start_idx:end_idx]
      
      batch_ts = data_tensors$time_series[batch_idx, ]
      batch_params = data_tensors$params[batch_idx, ]
      batch_scenarios = data_tensors$scenarios[batch_idx, ]
      
      output = model(batch_ts, batch_params, batch_scenarios)
      losses = vae_loss(output$reconstruction, batch_ts, 
                         output$mu, output$logvar, kl_weight)
      
      optimizer$zero_grad()
      losses$total$backward()
      optimizer$step()
      
      train_losses$total = train_losses$total + as.numeric(losses$total)
      train_losses$recon = train_losses$recon + as.numeric(losses$reconstruction)
      train_losses$kl = train_losses$kl + as.numeric(losses$kl)
    }
    
    # Validation phase
    model$eval()
    with_no_grad({
      val_output = model(val_ts, val_params, val_scenarios)
      val_losses = vae_loss(val_output$reconstruction, val_ts,
                             val_output$mu, val_output$logvar, kl_weight)
    })
    
    # Record history
    history$train_loss[epoch] = train_losses$total / n_batches
    history$val_loss[epoch] = as.numeric(val_losses$total)
    history$train_recon[epoch] = train_losses$recon / n_batches
    history$val_recon[epoch] = as.numeric(val_losses$reconstruction)
    history$train_kl[epoch] = train_losses$kl / n_batches
    history$val_kl[epoch] = as.numeric(val_losses$kl)
    
    # Early stopping
    if(history$val_loss[epoch] < best_val_loss) {
      best_val_loss = history$val_loss[epoch]
      patience_counter = 0
      # Save best model
      best_model_state = model$state_dict()
    } else {
      patience_counter = patience_counter + 1
    }
    
    if(patience_counter >= early_stopping_patience) {
      cat(sprintf("Early stopping at epoch %d\n", epoch))
      # Restore best model
      model$load_state_dict(best_model_state)
      break
    }
    
    if(epoch %% 10 == 0) {
      cat(sprintf("Epoch %d/%d - Train Loss: %.4f, Val Loss: %.4f\n",
                  epoch, epochs, history$train_loss[epoch], history$val_loss[epoch]))
    }
  }
  
  history
}

train_vae = function(model, data_tensors, epochs = 100, batch_size = 32, 
                     learning_rate = 0.001, kl_weight = 0.001) {
  
  optimizer = optim_adam(model$parameters, lr = learning_rate)
  
  n_samples = data_tensors$time_series$shape[1]
  n_batches = ceiling(n_samples / batch_size)
  
  history = list(
    total_loss = numeric(epochs),
    recon_loss = numeric(epochs),
    kl_loss = numeric(epochs)
  )
  
  for(epoch in 1:epochs) {
    model$train()
    epoch_losses = list(total = 0, recon = 0, kl = 0)
    
    # Shuffle indices
    indices = sample(n_samples)
    
    for(batch in 1:n_batches) {
      # Get batch indices
      start_idx = (batch - 1) * batch_size + 1
      end_idx = min(batch * batch_size, n_samples)
      batch_indices = indices[start_idx:end_idx]
      
      # Get batch data
      batch_ts = data_tensors$time_series[batch_indices, ]
      batch_params = data_tensors$params[batch_indices, ]
      batch_scenarios = data_tensors$scenarios[batch_indices, ]
      
      # Forward pass
      output = model(batch_ts, batch_params, batch_scenarios)
      
      # Calculate loss
      losses = vae_loss(output$reconstruction, batch_ts, 
                        output$mu, output$logvar, kl_weight)
      
      # Backward pass
      optimizer$zero_grad()
      losses$total$backward()
      optimizer$step()
      
      # Record losses
      epoch_losses$total = epoch_losses$total + as.numeric(losses$total)
      epoch_losses$recon = epoch_losses$recon + as.numeric(losses$reconstruction)
      epoch_losses$kl = epoch_losses$kl + as.numeric(losses$kl)
    }
    
    # Average losses
    history$total_loss[epoch] = epoch_losses$total / n_batches
    history$recon_loss[epoch] = epoch_losses$recon / n_batches
    history$kl_loss[epoch] = epoch_losses$kl / n_batches
    
    if(epoch %% 10 == 0) {
      cat(sprintf("Epoch %d/%d - Loss: %.4f (Recon: %.4f, KL: %.4f)\n",
                  epoch, epochs, history$total_loss[epoch],
                  history$recon_loss[epoch], history$kl_loss[epoch]))
    }
  }
  
  history
}

# ======================================================
# 4. ANALYSIS FUNCTIONS: UPDATED ANALYSIS FUNCTIONS FOR REPLICATES
# ======================================================

extract_latent_representations = function(model, data_tensors) {
  model$eval()
  
  with_no_grad({
    output = model(data_tensors$time_series, 
                   data_tensors$params, 
                   data_tensors$scenarios)
    
    latent_df = data.frame(
      as.matrix(output$z$cpu())
    )
    colnames(latent_df) = paste0("z", 1:ncol(latent_df))
    
    # Add metadata INCLUDING replicate_id
    cbind(latent_df, data_tensors$metadata)
  })
}

# UPDATED: Analyze Treg effects accounting for replicates
analyze_treg_effects_with_replicates = function(latent_df) {
  results = list()
  
  unique_params = unique(latent_df$param_id)
  
  for(param in unique_params) {
    param_data = latent_df %>% filter(param_id == param)
    
    for(sterile_cond in c(0, 1)) {
      # Get all replicates for each condition
      baseline_reps = param_data %>%
        filter(sterile == sterile_cond, allow_tregs == 0)
      
      treg_normal_reps = param_data %>%
        filter(sterile == sterile_cond, allow_tregs == 1, randomize_tregs == 0)
      
      treg_random_reps = param_data %>%
        filter(sterile == sterile_cond, allow_tregs == 1, randomize_tregs == 1)
      
      z_cols = grep("^z", colnames(baseline_reps))
      
      # Compare Tregs vs No Tregs (across all replicate pairs)
      if(nrow(baseline_reps) > 0 && nrow(treg_normal_reps) > 0) {
        distances_normal = c()
        for(i in 1:nrow(baseline_reps)) {
          for(j in 1:nrow(treg_normal_reps)) {
            dist = sqrt(sum((baseline_reps[i, z_cols] - treg_normal_reps[j, z_cols])^2))
            distances_normal = c(distances_normal, dist)
          }
        }
        
        results[[length(results) + 1]] = data.frame(
          param_id = param,
          sterile = sterile_cond,
          comparison = "treg_normal_vs_no_treg",
          mean_distance = mean(distances_normal),
          sd_distance = sd(distances_normal),
          median_distance = median(distances_normal),
          min_distance = min(distances_normal),
          max_distance = max(distances_normal),
          cv_distance = sd(distances_normal) / mean(distances_normal),  # Coefficient of variation
          n_comparisons = length(distances_normal)
        )
      }
      
      # Compare Randomized vs Normal Tregs
      if(nrow(treg_normal_reps) > 0 && nrow(treg_random_reps) > 0) {
        distances_random = c()
        for(i in 1:nrow(treg_normal_reps)) {
          for(j in 1:nrow(treg_random_reps)) {
            dist = sqrt(sum((treg_normal_reps[i, z_cols] - treg_random_reps[j, z_cols])^2))
            distances_random = c(distances_random, dist)
          }
        }
        
        results[[length(results) + 1]] = data.frame(
          param_id = param,
          sterile = sterile_cond,
          comparison = "treg_random_vs_treg_normal",
          mean_distance = mean(distances_random),
          sd_distance = sd(distances_random),
          median_distance = median(distances_random),
          min_distance = min(distances_random),
          max_distance = max(distances_random),
          cv_distance = sd(distances_random) / mean(distances_random),
          n_comparisons = length(distances_random)
        )
      }
    }
  }
  
  if(length(results) > 0) {
    bind_rows(results)
  } else {
    NULL
  }
}

# UPDATED: Parameter analysis using mean distances
identify_key_parameters = function(effects_df, parameters_df) {
  
  # Use mean_distance instead of latent_distance
  effects_with_params = effects_df %>%
    filter(comparison == "treg_normal_vs_no_treg") %>%
    left_join(parameters_df, by = c("param_id"="param_set_id"))
  
  # Calculate correlation between each parameter and mean Treg effect size
  param_correlations = list()
  
  for(param_name in parameter_names) {
    if(param_name %in% colnames(effects_with_params)) {
      cor_result = cor.test(
        effects_with_params[[param_name]], 
        effects_with_params$mean_distance,  # Changed from latent_distance
        method = "spearman"
      )
      
      param_correlations[[param_name]] = data.frame(
        parameter = param_name,
        correlation = cor_result$estimate,
        p_value = cor_result$p.value
      )
    }
  }
  
  correlations_df = bind_rows(param_correlations) %>%
    arrange(desc(abs(correlation)))
  
  # Also look at parameters that reduce variability (stabilizing effect)
  var_correlations = list()
  
  for(param_name in parameter_names) {
    if(param_name %in% colnames(effects_with_params)) {
      cor_result = cor.test(
        effects_with_params[[param_name]], 
        effects_with_params$cv_distance,  # Coefficient of variation
        method = "spearman"
      )
      
      var_correlations[[param_name]] = data.frame(
        parameter = param_name,
        variance_correlation = cor_result$estimate,
        variance_p_value = cor_result$p.value
      )
    }
  }
  
  variance_df = bind_rows(var_correlations) %>%
    arrange(desc(abs(variance_correlation)))
  
  list(
    correlations = correlations_df,
    variance_effects = variance_df,
    effects_summary = effects_with_params %>%
      summarize(
        mean_effect = mean(mean_distance),
        sd_effect = sd(mean_distance),
        mean_variability = mean(cv_distance)
      )
  )
}

# New visualization for replicate consistency
plot_replicate_consistency = function(effects_df) {
  effects_df %>%
    filter(comparison == "treg_normal_vs_no_treg") %>%
    ggplot(aes(x = mean_distance, y = cv_distance)) +
    geom_point(aes(color = factor(sterile)), alpha = 0.6) +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(
      title = "Treg Effect Consistency Across Replicates",
      x = "Mean Latent Distance (Effect Size)",
      y = "Coefficient of Variation (Replicate Consistency)",
      color = "Sterile"
    ) +
    annotate("text", x = Inf, y = 0.22, label = "High variability", 
             hjust = 1, color = "red")
}

# ======================================================
# 4. MAIN ANALYSIS SCRIPT - UPDATED
# ======================================================

### all in memory
run_treg_analysis = function(base_path, parameters_df, 
                             param_set_ids = NULL,
                             n_time_points = 500,
                             n_replicates = 10,  # ADD THIS
                             epochs = 100,
                             latent_dim = 32) {
  
  cat("Loading simulation data with all replicates...\n")
  
  if(is.null(param_set_ids)) {
    param_set_ids = unique(parameters_df$param_set_id)
  }
  
  # Load data WITH REPLICATES
  data_list = load_simulation_data(
    param_set_ids = param_set_ids[1:min(100, length(param_set_ids))],
    base_path = base_path,
    parameters_df = parameters_df,
    n_time_points = n_time_points,
    n_replicates = n_replicates  # ADD THIS
  )
  
  if(length(data_list) == 0) {
    stop("No data loaded!")
  }
  
  cat(sprintf("Loaded %d samples (including all replicates)\n", length(data_list)))
  cat(sprintf("Expected: %d params × 6 scenarios × %d replicates = %d\n", 
              length(param_set_ids[1:min(100, length(param_set_ids))]), 
              n_replicates, 
              length(param_set_ids[1:min(100, length(param_set_ids))]) * 6 * n_replicates))
  
  # Prepare tensors
  cat("Preparing tensors...\n")
  data_tensors = prepare_tensors(data_list, normalize = TRUE)
  
  # Create model
  cat("Creating VAE model...\n")
  model = VAE(
    time_steps = n_time_points,
    n_params = 24,
    n_scenarios = 6,
    latent_dim = latent_dim,
    hidden_dim = 256
  )
  
  # Train model (now with more samples due to replicates)
  cat("Training VAE (this will take longer with replicates)...\n")
  # history = train_vae(
  history = train_vae_with_validation(
    model = model,
    data_tensors = data_tensors,
    epochs = epochs,
    batch_size = 64,  # Might increase batch size
    learning_rate = 0.001,
    kl_weight = 0.001
  )
  
  # Extract latent representations
  cat("Extracting latent representations...\n")
  latent_df = extract_latent_representations(model, data_tensors)
  
  # Analyze Treg effects WITH REPLICATES
  cat("Analyzing Treg effects across replicates...\n")
  effects_df = analyze_treg_effects_with_replicates(latent_df)
  
  # Find parameters where Tregs are most beneficial
  if(!is.null(effects_df)) {
    top_beneficial = effects_df %>%
      filter(comparison == "treg_normal_vs_no_treg") %>%
      arrange(desc(mean_distance)) %>%
      head(20)
    
    cat("\nTop parameter sets where Tregs have largest effect:\n")
    print(top_beneficial %>% select(param_id, sterile, mean_distance, cv_distance))
    
    # Check consistency
    cat("\nReplicate consistency summary:\n")
    cat(sprintf("Mean CV across all parameters: %.3f\n", mean(effects_df$cv_distance)))
    cat(sprintf("Parameters with high variability (CV > 0.3): %d\n", 
                sum(effects_df$cv_distance > 0.3)))
  }
  
  # Return results
  list(
    model = model,
    history = history,
    latent_df = latent_df,
    effects_df = effects_df,
    data_tensors = data_tensors
  )
}

# Modified to process in chunks - batch loading 
run_treg_analysis_batched = function(base_path, parameters_df, 
                                      param_set_ids = NULL,
                                      batch_size = 10,  # Process 10 param sets at a time
                                      n_time_points = 500,
                                      n_replicates = 10,
                                      epochs = 100,
                                      latent_dim = 32) {
  
  if(is.null(param_set_ids)) {
    param_set_ids = unique(parameters_df$param_set_id)
  }
  
  # Process in batches
  n_batches = ceiling(length(param_set_ids) / batch_size)
  all_latent_dfs = list()
  
  # First, train the model on a subset
  cat("Training VAE on first batch...\n")
  first_batch_ids = param_set_ids[1:min(batch_size, length(param_set_ids))]
  
  # Load first batch
  data_list = load_simulation_data(
    param_set_ids = first_batch_ids,
    base_path = base_path,
    parameters_df = parameters_df,
    n_time_points = n_time_points,
    n_replicates = n_replicates
  )
  
  data_tensors = prepare_tensors(data_list, normalize = TRUE)
  
  # Create and train model
  model = VAE(time_steps = n_time_points, n_params = 24, 
               n_scenarios = 6, latent_dim = latent_dim)
  
  history = train_vae(model, data_tensors, epochs = epochs, 
                       batch_size = 64, learning_rate = 0.001)
  
  # Process all batches through trained model
  cat("Processing all parameter sets through trained VAE...\n")
  for(b in 1:n_batches) {
    start_idx = (b-1) * batch_size + 1
    end_idx = min(b * batch_size, length(param_set_ids))
    batch_ids = param_set_ids[start_idx:end_idx]
    
    cat(sprintf("Processing batch %d/%d (params %d-%d)\n", b, n_batches, start_idx, end_idx))
    
    # Load batch
    batch_data = load_simulation_data(
      param_set_ids = batch_ids,
      base_path = base_path,
      parameters_df = parameters_df,
      n_time_points = n_time_points,
      n_replicates = n_replicates
    )
    
    batch_tensors = prepare_tensors(batch_data, normalize = TRUE)
    
    # Extract latent representations
    latent_df = extract_latent_representations(model, batch_tensors)
    all_latent_dfs[[b]] = latent_df
    
    # Clear memory
    rm(batch_data, batch_tensors)
    gc()
  }
  
  # Combine all latent representations
  full_latent_df = bind_rows(all_latent_dfs)
  
  # Analyze effects
  effects_df = analyze_treg_effects_with_replicates(full_latent_df)
  
  list(
    model = model,
    history = history,
    latent_df = full_latent_df,
    effects_df = effects_df
  )
}
# ======================================================
# 5. EXECUTE ANALYSIS
# ======================================================

# Set paths
base_path = '/Users/burcutepekule/Desktop/tregs/mass_sim_results'
parameters_df = read_csv(
  '/Users/burcutepekule/Desktop/tregs/mass_sim_results/sampled_parameters.csv',
  show_col_types = FALSE
)

# Run analysis WITH REPLICATES
results = run_treg_analysis_batched(
  base_path = base_path,
  parameters_df = parameters_df,
  param_set_ids = 1:50,  # Start smaller since we have 10x more data
  batch_size    = 10,
  n_time_points = 500,
  n_replicates = 10,
  epochs = 50,
  latent_dim = 32
)

# Analyze parameter importance
if(!is.null(results$effects_df)) {
  param_analysis = identify_key_parameters(results$effects_df, parameters_df)
  
  cat("\n=== PARAMETER ANALYSIS ===\n")
  cat("\nTop 5 parameters correlated with mean Treg effectiveness:\n")
  print(param_analysis$correlations %>% head(5))
  
  cat("\nTop 5 parameters affecting replicate consistency:\n")
  print(param_analysis$variance_effects %>% head(5))
  
  # Plot replicate consistency
  p_consistency = plot_replicate_consistency(results$effects_df)
  print(p_consistency)
}

# Save results
torch_save(results$model, "treg_vae_model_with_replicates.pt")
saveRDS(results$latent_df, "latent_representations_with_replicates.rds")
saveRDS(results$effects_df, "treg_effects_with_replicates.rds")
cat("\nAnalysis complete!\n")
