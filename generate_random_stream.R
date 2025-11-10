#!/usr/bin/env Rscript
# Generate synchronized random number stream for reproducible scenario comparisons
rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)


total_draws_needed = 1e7

# Generate random numbers using R's RNG
for (int_in in 0:9999){
  set.seed(int_in)  # For reproducibility of the stream itself
  random_stream = runif(total_draws_needed)
  
  # Convert to data frame
  random_df = data.frame(random_number = random_stream)
  
  # Save to file
  output_file = paste0("./random_streams/random_numbers_seed_",int_in,".csv")
  write_csv(random_df, file = output_file)
  
  # cat(sprintf("\nRandom number stream saved to: %s\n", output_file))
  # cat(sprintf("File size: %.2f MB\n", file.info(output_file)$size / 1e6))
  # cat("\nUsage:\n")
  # cat("  - All simulations will draw from this stream sequentially\n")
  # cat("  - Reset the stream index before each new simulation\n")
  # cat("  - This ensures Tregs ON vs OFF use IDENTICAL random events\n")
}
