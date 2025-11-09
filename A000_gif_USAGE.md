# A000_gif_UPDATED.R Usage Guide

## Overview
Updated gif generation script that:
- Matches STANDALONE_SIM.R mechanics (all discrepancies fixed)
- Reads parameters from `balanced_lhs_parameters.csv`
- Generates simulation animations

## Key Updates from A000_gif.R

### 1. **DAMP Production**
- Now includes DAMPs from microbes touching epithelium (y=1)
- Uses logistic function for microbe-based DAMP contribution

### 2. **Event Ordering**
- Microbes move BEFORE DAMP calculation (not after)
- Ensures DAMPs reflect current microbe positions

### 3. **Epithelial Injury**
- Only pathogens cause injury (commensals removed)
- Uses `logistic_scaled_0_to_5_quantized()` instead of `log(count+1)`

### 4. **Treg Activation**
- Stochastic Beta-distributed antigen perception
- Uses `treg_discrimination_efficiency` parameter
- More realistic discrimination mechanism

## Command Line Usage

```bash
# Basic usage (parameter set 1, seed 3, infection, tregs ON, DAMPs-guided)
Rscript A000_gif_UPDATED.R

# Specify parameter set
Rscript A000_gif_UPDATED.R 5  # Use parameter set 5

# Full specification
Rscript A000_gif_UPDATED.R <param_row> <seed> <sterile> <allow_tregs> <randomize_tregs>

# Examples:
Rscript A000_gif_UPDATED.R 1 3 0 1 0  # Param 1, seed 3, infection, tregs ON, DAMPs-guided
Rscript A000_gif_UPDATED.R 10 42 1 1 0  # Param 10, seed 42, sterile, tregs ON, DAMPs-guided
Rscript A000_gif_UPDATED.R 5 7 0 0 1  # Param 5, seed 7, infection, tregs OFF, random movement
```

## Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `param_row` | 1 | Which row from balanced_lhs_parameters.csv to use |
| `seed` | 3 | Random seed for reproducibility |
| `sterile` | 0 | 0 = infection (pathogens), 1 = sterile injury |
| `allow_tregs` | 1 | 0 = tregs disabled, 1 = tregs active |
| `randomize_tregs` | 0 | 0 = follow DAMPs gradient, 1 = random movement |

## Output Files

### Generated Files
```
frames/
  frame_seed_<seed>_STERILE_<sterile>_TREGS_<allow_tregs>_trnd_<randomize_tregs>_paramset_<param_row>_<t>.png
  simulation_sterile<sterile>_tregs_<allow_tregs>_trnd_<randomize_tregs>_paramset_<param_row>.mp4
```

### Example Output Names
- `frames/frame_seed_3_STERILE_0_TREGS_1_trnd_0_paramset_1_100.png` - Frame 100
- `frames/simulation_sterile0_tregs_1_trnd_0_paramset_1.mp4` - Final video

## Parameters from CSV

The script reads these parameters from `balanced_lhs_parameters.csv`:

**Signal Thresholds:**
- `th_ROS_microbe` - ROS threshold for killing microbes
- `th_ROS_epith_recover` - ROS threshold for epithelial injury
- `activation_threshold_DAMPs` - DAMP threshold for phagocyte activation
- `activation_threshold_SAMPs` - SAMP threshold for phagocyte activation

**Diffusion & Decay:**
- `diffusion_speed_DAMPs`, `diffusion_speed_SAMPs`, `diffusion_speed_ROS`
- `DAMPs_decay`, `SAMPs_decay`, `ros_decay`

**Signal Production:**
- `add_DAMPs`, `add_SAMPs`, `add_ROS`

**Phagocyte Activity:**
- `activity_engulf_M0_baseline`, `activity_engulf_M1_baseline`, `activity_engulf_M2_baseline`
- `activity_ROS_M1_baseline`

**Leak Rates:**
- `rate_leak_pathogen_injury` - Pathogen influx rate at injury sites
- `rate_leak_commensal_injury` - Commensal influx at injury sites
- `rate_leak_commensal_baseline` - Baseline commensal influx

**Treg Parameters:**
- `treg_discrimination_efficiency` - How well Tregs discriminate antigens (0-1)
- `rat_com_pat_threshold` - Commensal ratio threshold for Treg activation

**Other:**
- `epith_recovery_chance` - Probability of epithelial recovery per timestep
- `active_age_limit` - Max age before cells can change phenotype

## Batch Processing Example

```bash
#!/bin/bash
# Run multiple parameter sets

for param_set in {1..10}; do
  for seed in 1 2 3; do
    echo "Running param_set=$param_set, seed=$seed"
    Rscript A000_gif_UPDATED.R $param_set $seed 0 1 0
  done
done
```

## Dependencies

Required R packages:
- `dplyr`
- `tidyr`
- `cowplot`
- `ggplot2`
- `gridExtra`
- `grid`
- `av` (for video encoding)

Required custom scripts:
- `MISC/FAST_FUNCTIONS.R`
- `MISC/PLOT_FUNCTIONS.R`
- `CONVERT_TO_DATAFRAME.R`

## Notes

- Video encoding requires `av` package and ffmpeg
- Frame generation can be slow (600 dpi images)
- Each simulation runs for 500 timesteps
- Grid size is fixed at 25x25
- Agent counts: 35% of grid for both phagocytes and Tregs
