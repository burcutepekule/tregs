#!/usr/bin/env Rscript
# Generate synchronized random number stream for reproducible scenario comparisons
# This ensures Tregs ON vs OFF comparisons use identical random events

library(here)

# Set parameters
n_simulations <- 100  # How many full simulations to support
t_max <- 500          # Timesteps per simulation
safety_factor <- 10   # Multiply by this to ensure enough random numbers

# Estimate random draws per timestep (conservative upper bound)
grid_size <- 25
n_agents <- round(grid_size^2 * 0.35) * 2  # phagocytes + tregs
n_microbes_max <- 500  # conservative max

draws_per_timestep <- (
  n_agents * 2 +           # agent movement (2 draws per agent for dx, dy)
  n_microbes_max * 2 +     # microbe movement
  n_microbes_max * 2 +     # ROS killing checks
  n_agents * 10 +          # engulfment, activation, etc
  grid_size * 2 +          # epithelial injury/recovery
  100                      # buffer for misc operations
)

total_draws_needed <- n_simulations * t_max * draws_per_timestep * safety_factor

cat(sprintf("Generating %d random numbers...\n", total_draws_needed))
cat(sprintf("  This supports ~%d simulations of %d timesteps\n", n_simulations, t_max))

# Generate random numbers using R's high-quality RNG
set.seed(42)  # For reproducibility of the stream itself
random_stream <- runif(total_draws_needed)

# Save to file
output_file <- here("random_numbers_seed_42.txt")
write.table(random_stream,
            file = output_file,
            row.names = FALSE,
            col.names = FALSE)

cat(sprintf("\nRandom number stream saved to: %s\n", output_file))
cat(sprintf("File size: %.2f MB\n", file.info(output_file)$size / 1e6))
cat("\nUsage:\n")
cat("  - All simulations will draw from this stream sequentially\n")
cat("  - Reset the stream index before each new simulation\n")
cat("  - This ensures Tregs ON vs OFF use IDENTICAL random events\n")
