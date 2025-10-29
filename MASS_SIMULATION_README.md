# Mass Simulation Framework

Parallel simulation framework for running thousands of ABM simulations with sampled parameters.

## Quick Start

### 1. Test the Framework (2 minutes)

```bash
python test_mass_sim.py
```

This runs a small test with:
- 2 parameter sets
- 2 replicates
- 2 cores
- Should complete in ~1-2 minutes

### 2. Run Full Mass Simulation

```bash
# Run 100 parameter sets with 100 replicates each (10 cores)
python mass_simulation.py --n_param_sets 100 --n_replicates 100 --n_cores 10
```

**Expected time:**
- ~3 seconds per simulation
- 100 param sets × 6 scenarios × 100 replicates = 60,000 simulations
- 60,000 sims × 3 sec / 10 cores = **~5 hours**

### 3. Analyze Results

```bash
python analyze_mass_sim.py --results_dir mass_sim_results
```

This creates:
- `analysis_full_with_params.csv` - All replicates with parameters
- `analysis_aggregated.csv` - Mean ± std per parameter-scenario combo

---

## Parameter Sampling

Parameters are sampled uniformly from these ranges:

### Thresholds
- `th_ROS_microbe`: [0, 0.5] - ROS level that kills bacteria
- `th_ROS_epith_recover`: [th_ROS_microbe, 1] - ROS level that injures epithelium
- `epith_recovery_chance`: [0, 1] - Probability of epithelial cell recovery
- `rat_com_pat_threshold`: [0.5, 1] - Threshold for Treg activation

### Diffusion (MUST be < 0.125 for stability!)
- `diffusion_speed_DAMPs`: [0, 0.12]
- `diffusion_speed_SAMPs`: [0, 0.12]
- `diffusion_speed_ROS`: [0, 0.12]

### Signal Production
- `add_ROS`: [0, 1]
- `add_DAMPs`: [0, 1]
- `add_SAMPs`: [0, 1]

### Decay Rates
- `DAMPs_decay`: [0, 1]
- `SAMPs_decay`: [0, 1]
- `ros_decay`: [0, 1]

### Activation Thresholds
- `activation_threshold_DAMPs`: [0, 1]
- `activation_threshold_SAMPs`: [0, 1]

### Activity Levels (max 50% of capacity)
- `activity_engulf_M0_baseline`: [0, 0.5]
- `activity_engulf_M1_baseline`: [0, 0.5]
- `activity_engulf_M2_baseline`: [0, 0.5]
- `activity_ROS_M1_baseline`: [0, 0.5]

### Leakage Rates
- `rate_leak_commensal_injury`: [0.5, 1]
- `rate_leak_pathogen_injury`: [0.5, 1]
- `rate_leak_commensal_baseline`: [0, 0.25]

### Other
- `active_age_limit`: {3, 4, ..., 30} - How long cells stay activated
- `treg_discrimination_efficiency`: [0, 1] - How well Tregs distinguish commensals from pathogens

---

## Scenarios

6 valid scenarios are generated:

| ID | Sterile | Allow Tregs | Suppress Cognate | Randomize Tregs |
|----|---------|-------------|------------------|-----------------|
| 0  | No      | No          | No               | No              |
| 1  | No      | Yes         | No               | No              |
| 2  | No      | Yes         | No               | Yes             |
| 3  | Yes     | No          | No               | No              |
| 4  | Yes     | Yes         | No               | No              |
| 5  | Yes     | Yes         | No               | Yes             |

**Note:** Scenarios with `allow_tregs=No` and `randomize=Yes` are filtered out (doesn't make sense).

---

## File Structure

After running mass simulations, you'll have:

```
mass_sim_results/
├── sampled_parameters.csv          # Parameter sets (one row per param_set_id)
├── scenarios.csv                   # Scenario definitions (one row per scenario_id)
├── all_simulation_results.csv      # Raw results (all timepoints, all replicates)
├── analysis_full_with_params.csv   # Summary + parameters merged
└── analysis_aggregated.csv         # Mean ± std per parameter-scenario combo
```

---

## Output Columns

### all_simulation_results.csv
Same structure as `simulation_results.csv`, plus:
- `param_set_id`: Which parameter set (0 to n_param_sets-1)
- `scenario_id`: Which scenario (0 to 5)
- `replicate_id`: Which replicate (0 to 99)
- `seed`: Random seed used
- `t`: Timestep (0 to 500)
- All agent counts (epithelial, bacteria, phagocytes, tregs, etc.)

### analysis_full_with_params.csv
One row per simulation (param_set × scenario × replicate):
- All identifiers (param_set_id, scenario_id, replicate_id)
- All sampled parameters
- All scenario settings
- Summary metrics:
  - `pathogen_clearance_time`: When pathogens dropped to zero
  - `pathogen_peak`: Maximum pathogen count
  - `commensal_final`: Final commensal count
  - `epithelial_healthy_final`: Final healthy epithelial cells
  - `M1_peak`, `M2_peak`: Peak phagocyte phenotypes
  - `treg_active_mean`: Average active Treg count
  - `health_score`: Overall health metric (0=bad, 1=good)

### analysis_aggregated.csv
One row per parameter-scenario combo:
- All parameters and scenario settings
- Mean and std for each metric across 100 replicates
- Example: `pathogen_clearance_time_mean`, `pathogen_clearance_time_std`

---

## Advanced Options

### Custom Parameter Ranges

Edit `mass_simulation.py`, function `sample_parameters()`:

```python
# Example: Narrow down Treg efficiency
treg_discrimination_efficiency = np.random.uniform(0.7, 1.0, n_sets)  # Only high values
```

### Custom Scenarios

Edit `mass_simulation.py`, function `generate_scenarios()`:

```python
# Example: Add allow_suppress=True scenarios
for allow_suppress in [0, 1]:  # Instead of just [0]
    scenarios.append({...})
```

### Checkpoint/Resume (for very long runs)

The current implementation runs all simulations in one go. For multi-day runs, consider:

1. Split into batches:
```bash
python mass_simulation.py --n_param_sets 100 --output_dir batch_1
python mass_simulation.py --n_param_sets 100 --output_dir batch_2 --base_seed 1000
```

2. Combine results later:
```python
import pandas as pd
df1 = pd.read_csv('batch_1/all_simulation_results.csv')
df2 = pd.read_csv('batch_2/all_simulation_results.csv')
combined = pd.concat([df1, df2], ignore_index=True)
combined.to_csv('combined_results.csv', index=False)
```

---

## Memory Management

Each simulation produces ~500 rows (time steps) × ~40 columns = ~20KB.

**Total memory for 60,000 simulations:**
- Results: 60,000 × 20KB = 1.2 GB
- Parameters: ~100 rows × 25 columns = negligible
- Summary: 60,000 × 30 columns = ~7 MB

**Should easily fit in RAM** on any modern machine.

If memory is a concern, modify `run_mass_simulations()` to save results in batches:
```python
# Save every 1000 simulations
if len(results_list) >= 1000:
    pd.concat(results_list).to_csv(f'batch_{batch_num}.csv')
    results_list = []
    batch_num += 1
```

---

## Troubleshooting

### "Too many open files" error
Increase file descriptor limit:
```bash
ulimit -n 4096
```

### Slow performance
- Check CPU usage: `htop` or Task Manager
- Reduce `n_cores` if system is overloaded
- Close other applications

### Out of memory
- Reduce `n_param_sets`
- Run in batches (see Advanced Options)
- Increase system RAM

### Results don't match expectations
- Check parameter ranges in `sample_parameters()`
- Verify scenario filtering in `generate_scenarios()`
- Look at `sampled_parameters.csv` to see what was sampled

---

## Next Steps

After analysis, you can:

1. **Import to R for visualization:**
```r
library(tidyverse)
results <- read_csv("mass_sim_results/analysis_aggregated.csv")
```

2. **Run sensitivity analysis:**
- Correlate parameters with health outcomes
- Identify key parameters driving system behavior

3. **Run focused follow-ups:**
- Identify interesting parameter regions
- Run more replicates in those regions

4. **Visualize parameter space:**
- Plot health_score vs parameters
- Identify parameter combinations that lead to good/bad outcomes
