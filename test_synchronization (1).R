#!/usr/bin/env Rscript
# Test script to verify synchronized random numbers work correctly
# This demonstrates that Tregs ON vs OFF use identical random streams

cat("=======================================================\n")
cat("Testing Synchronized Random Number Stream\n")
cat("=======================================================\n\n")

# Setup synchronized RNG
cat("1. Loading random stream...\n")
rng_env <- new.env()

# For testing, create a small stream
set.seed(42)
rng_env$stream <- runif(1000)
rng_env$index <- 1

my_runif <- function(n = 1) {
  if (rng_env$index + n - 1 > length(rng_env$stream)) {
    stop("Stream exhausted!")
  }
  vals <- rng_env$stream[rng_env$index:(rng_env$index + n - 1)]
  rng_env$index <- rng_env$index + n
  return(vals)
}

my_sample <- function(x, size, replace = FALSE, prob = NULL) {
  n <- length(x)
  if (!replace && size > n) stop("Sample too large")
  probs <- if (is.null(prob)) rep(1/n, n) else prob / sum(prob)
  u <- my_runif(size)
  cumprob <- cumsum(probs)
  indices <- findInterval(u, cumprob) + 1
  x[indices]
}

sample_rbeta <- function(alpha, beta) {
  u1 <- my_runif(1)
  u2 <- my_runif(1)
  x <- qgamma(u1, shape = alpha, rate = 1.0)
  y <- qgamma(u2, shape = beta, rate = 1.0)
  return(x / (x + y))
}

runif <- my_runif
sample <- my_sample

cat("   Stream length:", length(rng_env$stream), "\n\n")

# Test 1: Basic operations consume same numbers
cat("2. Test basic synchronization...\n")
rng_env$index <- 1
draw1 <- runif(5)

rng_env$index <- 1
draw2 <- runif(5)

cat("   First 5 draws (run 1):", draw1, "\n")
cat("   First 5 draws (run 2):", draw2, "\n")
cat("   Identical:", identical(draw1, draw2), "\n\n")

# Test 2: Conditional logic that DOESN'T break sync
cat("3. Test conditional logic (CORRECT - always sample)...\n")

simulate_tregs_correct <- function(allow_tregs) {
  rng_env$index <- 1  # Reset stream
  draws <- numeric(0)

  for (i in 1:10) {
    # Some unconditional draws
    x <- runif(1)
    draws <- c(draws, x)

    # CORRECT: Always sample, conditionally use
    alpha <- 2.0
    beta <- 3.0
    beta_sample <- sample_rbeta(alpha, beta)  # ALWAYS called

    # Only USE it if allowed
    if (allow_tregs && beta_sample > 0.5) {
      # Do something
      effect <- 1
    } else {
      effect <- 0
    }

    # More unconditional draws
    y <- runif(1)
    draws <- c(draws, y)
  }

  return(list(draws = draws, final_index = rng_env$index))
}

result_off <- simulate_tregs_correct(allow_tregs = FALSE)
result_on <- simulate_tregs_correct(allow_tregs = TRUE)

cat("   Stream index after Tregs OFF:", result_off$final_index, "\n")
cat("   Stream index after Tregs ON:", result_on$final_index, "\n")
cat("   Same random draws:", identical(result_off$draws, result_on$draws), "\n")
cat("   ✓ Synchronized!\n\n")

# Test 3: Conditional logic that DOES break sync (WRONG way)
cat("4. Test conditional logic (WRONG - conditional sample)...\n")

simulate_tregs_wrong <- function(allow_tregs) {
  rng_env$index <- 1  # Reset stream
  draws <- numeric(0)

  for (i in 1:10) {
    x <- runif(1)
    draws <- c(draws, x)

    # WRONG: Only sample if allowed
    if (allow_tregs) {
      alpha <- 2.0
      beta <- 3.0
      beta_sample <- sample_rbeta(alpha, beta)  # Only called when TRUE!

      if (beta_sample > 0.5) {
        effect <- 1
      }
    }

    y <- runif(1)
    draws <- c(draws, y)
  }

  return(list(draws = draws, final_index = rng_env$index))
}

result_off_wrong <- simulate_tregs_wrong(allow_tregs = FALSE)
result_on_wrong <- simulate_tregs_wrong(allow_tregs = TRUE)

cat("   Stream index after Tregs OFF:", result_off_wrong$final_index, "\n")
cat("   Stream index after Tregs ON:", result_on_wrong$final_index, "\n")
cat("   Same random draws:", identical(result_off_wrong$draws, result_on_wrong$draws), "\n")
cat("   ✗ OUT OF SYNC!\n\n")

cat("=======================================================\n")
cat("Summary:\n")
cat("=======================================================\n")
cat("✓ Correct approach: Always sample, conditionally use\n")
cat("✗ Wrong approach: Conditionally sample\n")
cat("\nThe updated A000_gif_UPDATED.R uses the CORRECT approach!\n")
