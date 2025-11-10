# Logistic function
logistic_function <- function(x, k = k_in, x0 = x0_in) {
  return(1 / (1 + exp(-k * (x - x0))))
}

logistic_scaled_0_to_5_quantized <- function(x,  k = k_in, x0 = x0_in) {
  # return(round(4*plogis(x, location = x0, scale = 1 / k))+1)
  return(round(5*plogis(x, location = x0, scale = 1 / k)))
}

# k_in = 0.044
# x0_in= 50
# plot(seq(0,200,1),round(log(1+seq(0,200,1))),col='red')
# points(seq(0,200,1),logistic_scaled_0_to_5_quantized(seq(0,200,1),k=k_in,x0=x0_in), col='blue')

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

# Optimized version of iszero_coordinates
iszero_coordinates <- function(x) {
  # Initialize with default sampling: -1, 0, 1
  y <- sample(c(-1, 0, 1), length(x), replace = TRUE)
  # Replace where x == 0 with sample from -1 or 1
  zero_idx <- which(x == 0)
  y[zero_idx] <- sample(c(-1, 1), length(zero_idx), replace = TRUE)
  
  return(y)
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
shift_insert <- function(current_registry, new_elements_vector) {
  combined_registry <- c(new_elements_vector, current_registry)
  target_length <- length(current_registry)
  result_registry <- combined_registry[1:target_length]
  
  return(result_registry)
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

sample_rbeta <- function(alpha, beta) {
  x <- rgamma(1, shape = alpha, rate = 1.0)
  y <- rgamma(1, shape = beta, rate = 1.0)
  return(x / (x + y))
}


