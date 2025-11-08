# Analysis: Why You Get ~4% Instead of ~25-30% "Worse" Outcomes

## Problem Summary

You expected ~25-30% "worse" outcomes based on your decision tree balancing (line 300 in A006_decision_tree_a_bias_class.R), but you're only getting ~4% when you run simulations with the pre-sampled parameters.

## Root Cause

### The Flawed Sampling Strategy (Lines 378-410)

Your R script does the following for each decision tree node:

```r
# 1. Get parameter bounds (5th-95th percentile of all data in this node)
bounds = get_node_bounds(df_clustering, node_id, param_names)

# 2. Generate samples for ALL categories from THE SAME bounds
samples_better = generate_lhs_samples(bounds, n_better, param_names)  # ← same bounds
samples_worse = generate_lhs_samples(bounds, n_worse, param_names)   # ← same bounds
samples_drift = generate_lhs_samples(bounds, n_drift, param_names)   # ← same bounds

# 3. Label them by category
samples_better$target_category = "better"  # ← just a label, doesn't affect sampling!
samples_worse$target_category = "worse"    # ← just a label!
```

**The fundamental flaw:** All three categories are sampled **uniformly** (via LHS) from the **same parameter bounds**. There's no mechanism to ensure that samples labeled "worse" actually come from parameter regions that produce worse outcomes.

### Why This Fails

1. **Uniform sampling within bounds**: LHS ensures good parameter space coverage but samples uniformly. If "worse" outcomes came from a specific sub-region within a node's parameter space, uniform sampling will dilute them.

2. **Quantile-based bounds exclude extremes**: Using 5th-95th percentiles removes 10% of original data, potentially excluding extreme parameter values that strongly correlate with specific outcomes.

3. **Labeling ≠ Targeting**: The script labels samples as "worse" based on the expectation that they'll follow the node's probability distribution, but there's no actual targeting mechanism.

### The Math

Your expected proportion calculation (line 300):
```r
worse_rat = sum(node_probabilities$worse * samples_per_node) / total_samples
```

This assumes: **P(worse | new_sample_from_node_X) = P(worse | original_data_in_node_X)**

But this breaks down because:
- Original data: Actual simulation results from specific parameter combinations
- New LHS samples: Uniform coverage of parameter bounds without outcome targeting
- **P(worse | uniform_LHS_in_bounds) ≠ P(worse | original_data)**

## Solutions

I've created two fixed versions:

### Option 1: Direct Resampling (Recommended for Accuracy)

**File:** `A006_decision_tree_FIXED_direct_resample.R`

**Strategy:**
- Directly resample from the original parameter sets that produced each outcome
- Sample with replacement if needed to achieve target counts
- Guarantees outcome proportions match expectations

**Pros:**
- Simple and reliable
- Guaranteed to produce the expected outcome proportions
- No loss of information

**Cons:**
- Limited exploration of new parameter space
- May oversample some parameter sets if using replacement

### Option 2: Perturbed Sampling (Recommended for Exploration)

**File:** `A006_decision_tree_FIXED_perturbed.R`

**Strategy:**
- Start with original parameter sets that produced each outcome
- Add small random perturbations (5% of parameter range)
- Explore nearby parameter space while staying close to known outcomes

**Pros:**
- Explores new parameter space
- Better coverage than direct resampling
- Still maintains strong connection to known outcomes

**Cons:**
- Small chance that perturbed parameters produce different outcomes
- Slightly less guaranteed outcome proportions

## Key Changes in Fixed Versions

1. **Separate by outcome FIRST**, then sample:
   ```r
   params_better = df_clustering %>% filter(dominant_outcome == "better")
   params_worse = df_clustering %>% filter(dominant_outcome == "worse")
   params_drift = df_clustering %>% filter(dominant_outcome == "drift")
   ```

2. **Keep target_category for validation**:
   ```r
   # OLD: balanced_lhs_dataset[c('param_set_id', param_names)]
   # NEW: Keep target_category column in saved CSV
   write.csv(balanced_lhs_dataset, "balanced_lhs_parameters_FIXED.csv")
   ```

3. **Sample from outcome-specific sets**:
   - Not from node bounds
   - Not with LHS within bounds
   - Directly from parameter sets that produced each outcome

## How to Use

### Step 1: Generate Balanced Parameters

Choose one approach:

```bash
# Option 1: Direct resampling (most reliable)
Rscript A006_decision_tree_FIXED_direct_resample.R

# Option 2: Perturbed sampling (better exploration)
Rscript A006_decision_tree_FIXED_perturbed.R
```

This will create:
- `balanced_lhs_parameters_DIRECT.csv` or
- `balanced_lhs_parameters_PERTURBED.csv`

**Important:** These files now INCLUDE the `target_category` column for validation!

### Step 2: Run Simulations

```bash
python mass_simulation_LHS_presampled.py \
  --csv_file balanced_lhs_parameters_DIRECT.csv \
  --n_cores 10 \
  --output_dir mass_sim_results_presampled_FIXED
```

### Step 3: Analyze and Validate

After running simulations and computing outcomes (A008, A009), you can validate:

```r
# Load results
results = readRDS('all_comparison_results_0_pre.rds')
params = read_csv('balanced_lhs_parameters_DIRECT.csv')

# Join and compare expected vs actual
validation = results %>%
  left_join(params %>% select(param_set_id, target_category), by='param_set_id') %>%
  mutate(actual_outcome = ...) # your outcome classification

# Confusion matrix
table(validation$target_category, validation$actual_outcome)
```

## Expected Results

With the fixed approach:

- **Direct resampling:** Should see ~33% better, ~33% worse, ~33% drift
- **Perturbed sampling:** Should see close to ~33% each, with small variations due to perturbations

## Additional Notes

### Why the Original Approach Seemed Reasonable

The decision tree provides node-level probabilities, which is valuable information. The mistake was assuming that:
1. These probabilities would transfer to new LHS samples within node bounds
2. Uniform sampling within bounds preserves outcome distributions

This works for **parameter bounds** but fails for **outcome targeting**.

### Better Alternatives (Future Work)

For more sophisticated balancing, consider:

1. **Conditional VAE:** Train a generative model on each outcome class separately
2. **SMOTE-like approaches:** Synthesize minority class samples
3. **Refined decision tree rules:** Use the tree to identify specific parameter combinations, not just bounds
4. **Stratified sampling:** Sample within finer subdivisions of each node based on parameter interactions

## Summary

**The issue:** Sampling uniformly within node bounds ≠ sampling with outcome targeting

**The fix:** Sample from parameter sets that actually produced each outcome

**The result:** Balanced outcomes that match your expectations
