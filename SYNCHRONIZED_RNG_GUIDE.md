# Synchronized Random Numbers for Fair Scenario Comparisons

## The Problem

When comparing scenarios (e.g., Tregs ON vs OFF), you want to know:
- **Is the difference due to the Treg mechanism?**
- **Or just random variation between runs?**

With standard `set.seed()`, each run gets different random numbers, making it impossible to isolate the biological effect from stochastic noise.

## The Solution: Common Random Numbers (CRN)

Use the **same sequence** of random numbers for all scenarios:

```
Scenario A (Tregs OFF):  R1, R2, R3, R4, R5, ...
Scenario B (Tregs ON):   R1, R2, R3, R4, R5, ...  ← SAME sequence!
```

This way:
- Initial positions are identical
- Microbe movements are identical
- All stochastic events are identical
- **Only the Treg mechanism differs**

## Setup (One-Time)

### Step 1: Generate Random Stream

```bash
Rscript generate_random_stream.R
```

This creates `random_numbers_seed_42.txt` (~millions of pre-generated random numbers)

### Step 2: Verify File Created

```bash
ls -lh random_numbers_seed_42.txt
# Should show ~50-100 MB file
```

## Usage: Fair Comparisons

### Example 1: Tregs ON vs OFF

```bash
# Run BOTH with seed_in=1 (same random stream)
Rscript A000_gif_UPDATED.R 1 1 0 0 0  # Param 1, seed 1, infection, Tregs OFF
Rscript A000_gif_UPDATED.R 1 1 0 1 0  # Param 1, seed 1, infection, Tregs ON
#                          ↑ ↑ ↑ ↑ ↑
#                          │ │ │ │ └─ randomize_tregs
#                          │ │ │ └─── allow_tregs (ONLY DIFFERENCE!)
#                          │ │ └───── sterile
#                          │ └─────── seed (MUST BE SAME!)
#                          └───────── param_row
```

**Result:** Differences are **100% due to Treg mechanism**, not randomness!

### Example 2: Sterile vs Infection

```bash
# Compare sterile injury vs infection with Tregs ON
Rscript A000_gif_UPDATED.R 5 42 1 1 0  # Sterile injury
Rscript A000_gif_UPDATED.R 5 42 0 1 0  # Infection
#                             ↑↑       # SAME seed for both!
```

### Example 3: Multiple Replicates

For statistical power, run multiple seeds:

```bash
#!/bin/bash
# Compare Tregs ON vs OFF across 10 replicates

for seed in {1..10}; do
  echo "Replicate $seed"
  Rscript A000_gif_UPDATED.R 1 $seed 0 0 0  # Tregs OFF
  Rscript A000_gif_UPDATED.R 1 $seed 0 1 0  # Tregs ON
done

# Each seed gives a PAIRED comparison:
# - Seed 1: Tregs OFF vs ON (same randomness)
# - Seed 2: Tregs OFF vs ON (same randomness)
# - etc.
```

## How It Works

1. **Pre-generated stream**: File contains millions of random numbers
2. **Sequential drawing**: Each `runif()`, `sample()`, etc. pulls from the stream
3. **Synchronized runs**: Same seed → start at same position → identical draws
4. **Isolation**: Only mechanistic changes create differences

### Critical Implementation Detail: Always Consume Random Numbers

**THE GOLDEN RULE:** Always call random number generators, even when the result isn't used.

#### ❌ WRONG (breaks synchronization):
```r
if (allow_tregs_to_do_their_job) {
  rat_com_pat = sample_rbeta(alpha, beta)  # Only called when TRUE!
  if (rat_com_pat > threshold) {
    // activate tregs
  }
}
```

**Problem:**
- Tregs OFF: skips `sample_rbeta()`, uses draws [R1, R2, R3, ...]
- Tregs ON: calls `sample_rbeta()`, uses draws [R1, R2, **Rbeta**, R4, ...]
- Streams are now OUT OF SYNC!

#### ✅ CORRECT (maintains synchronization):
```r
// ALWAYS sample (both scenarios)
rat_com_pat = sample_rbeta(alpha, beta)

// ONLY use the result if allowed
if (allow_tregs_to_do_their_job && rat_com_pat > threshold) {
  // activate tregs
}
```

**Result:**
- Tregs OFF: calls `sample_rbeta()` but ignores result
- Tregs ON: calls `sample_rbeta()` and uses result
- Both use draws [R1, R2, **Rbeta**, R4, ...] → SYNCHRONIZED!

## Technical Details

### Custom RNG Functions

The script overrides R's built-in functions:
- `runif()` → draws from pre-generated stream
- `sample()` → uses stream via `runif()`
- `rbeta()` → uses stream for beta sampling (Treg discrimination)
- `rgamma()` → uses stream for gamma sampling

### Stream Format

```
0.12345
0.67890
0.45678
...
```

Simple text file, one number per line, values in [0,1].

### Stream Exhaustion

If you run out of random numbers:
```
Error: Random stream exhausted! Used 1234567 numbers, need 10 more.
```

**Solution:** Re-run `generate_random_stream.R` with higher `n_simulations`

## Verification: Are Comparisons Fair?

### Test 1: Identical Initial Conditions

```r
# Check initial microbe positions are identical
seed <- 1
# Run 1
# Run 2 (same seed)
# → Initial pathogen_coords should be IDENTICAL
```

### Test 2: Divergence Only from Mechanism

```r
# Tregs OFF vs ON (same seed)
# → Early timesteps: nearly identical
# → Later: diverge ONLY due to Treg-mediated SAMPs
```

### Test 3: Reproducibility

```r
# Run simulation twice with same seed
# → Should produce BITWISE IDENTICAL results
```

## Best Practices

### ✅ DO:
- Use **same seed** for direct comparisons (Tregs ON vs OFF)
- Use **different seeds** for independent replicates
- Generate **enough random numbers** (default supports 100+ full runs)
- Document which seed was used for each run

### ❌ DON'T:
- Compare runs with different seeds (not fair!)
- Mix synchronized and standard RNG in same analysis
- Forget to regenerate stream if you change parameters

## Disabling Synchronized RNG

If you want standard R randomness:

```bash
# Method 1: Command line argument
Rscript A000_gif_UPDATED.R 1 42 0 1 0 FALSE
#                                      ↑↑↑↑↑
#                                      use_synchronized_rng=FALSE

# Method 2: Edit script
# Set: use_synchronized_rng = FALSE (line 21)
```

## Statistical Analysis Example

```r
library(tidyverse)

# Load paired results
results <- tibble(
  seed = rep(1:10, each = 2),
  tregs = rep(c("OFF", "ON"), 10),
  final_pathogens = c(...)  # Load from your output
)

# Paired t-test (each seed is a matched pair)
results_wide <- results %>%
  pivot_wider(names_from = tregs, values_from = final_pathogens)

t.test(results_wide$OFF, results_wide$ON, paired = TRUE)

# This tests: Does Treg mechanism reduce pathogens?
# (controlling for random variation via pairing)
```

## Advanced: Multiple Scenarios

Compare 4 scenarios with full factorial design:

```bash
for seed in {1..5}; do
  for sterile in 0 1; do
    for tregs in 0 1; do
      Rscript A000_gif_UPDATED.R 1 $seed $sterile $tregs 0
    done
  done
done

# Creates 5 replicates × 4 scenarios = 20 runs
# Each replicate (seed) uses SAME random stream for all 4 scenarios
```

## Troubleshooting

### "Stream exhausted" error
**Solution:** Increase `n_simulations` in `generate_random_stream.R`

### Different results with same seed
**Possible causes:**
1. Different parameter values
2. Code version changed
3. Stream file changed
4. Synchronized RNG disabled

### Results still too variable
**Causes:**
- Comparing different seeds (not paired!)
- Need more replicates
- System is genuinely chaotic (biological reality!)

## References

This technique is called **Common Random Numbers (CRN)** and is standard practice in stochastic simulation:

- Law, A.M. (2015). Simulation Modeling and Analysis. McGraw-Hill.
- Nelson, B.L. (2013). Foundations and Methods of Stochastic Simulation. Springer.
