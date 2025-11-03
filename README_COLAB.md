# Running Mass Simulations on Google Colab

This guide will help you run your mass simulations on Google Colab with maximum CPU parallelization using your purchased compute units.

## 📁 Files for Google Colab

The following files are specifically designed for Google Colab (with `_G` suffix):

- **`mass_simulation_LHS_G.py`** - Main simulation runner with Colab optimizations
- **`main_simulation_G.py`** - Core simulation logic (Colab-friendly paths)
- **`simulation_utils_G.py`** - Utility functions with Numba JIT compilation
- **`colab_setup.ipynb`** - Interactive notebook for easy setup and execution

## 🚀 Quick Start Guide

### Step 1: Open Google Colab
1. Go to [Google Colab](https://colab.research.google.com/)
2. Sign in with your Google account
3. Click "File" → "Upload notebook"
4. Upload `colab_setup.ipynb`

### Step 2: Enable Compute Units
Since you've purchased compute units from Google Colab:
1. Click "Runtime" → "Change runtime type"
2. Select the best available hardware (with your compute units, you should have access to more CPUs)
3. Click "Save"

### Step 3: Follow the Notebook
The notebook has clear instructions for each step:
1. Check available resources (CPU cores, RAM)
2. Install required packages
3. Upload the 3 Python files (`*_G.py`)
4. Configure simulation parameters
5. Run the simulation
6. Download results

## ⚙️ Configuration Options

### Key Parameters

In the notebook, you can adjust:

```python
N_PARAM_SETS = 100      # Number of parameter sets (LHS sampling)
N_REPLICATES = 10       # Replicates per parameter-scenario combo
N_CORES = None          # CPU cores (None = use all available)
OUTPUT_DIR = 'mass_sim_results'
RANDOM_SEED = 42
COMBINE_RESULTS = False # Combine all results into one file?
```

### What Gets Generated

The simulation will create:
- `sampled_parameters.csv` - All sampled parameter sets
- `scenarios.csv` - All scenario combinations
- `simulation_results_param_set_N.csv` - Results for each parameter set

Total simulations = `N_PARAM_SETS × 6 scenarios × N_REPLICATES`

## 💡 Optimization Tips

### Maximize CPU Usage
- The scripts automatically detect and use all available CPU cores
- With compute units, you should have access to more powerful CPUs
- The simulation shows real-time progress and ETA

### Memory Considerations
- Results are saved incrementally (one file per parameter set)
- If you run out of memory, reduce `N_PARAM_SETS` or `N_REPLICATES`
- Set `COMBINE_RESULTS = False` to save memory

### Session Management
- Google Colab sessions can timeout after a few hours of inactivity
- If your session times out, simply re-run the notebook
- The simulation will automatically resume from where it left off (skips completed parameter sets)

## 📊 Understanding the Output

### File Structure
```
mass_sim_results/
├── sampled_parameters.csv              # Parameter sets (LHS sampled)
├── scenarios.csv                       # Scenario combinations
├── simulation_results_param_set_0.csv  # Results for param set 0
├── simulation_results_param_set_1.csv  # Results for param set 1
└── ...
```

### Result Columns
Each result file contains:
- Time series data (500 time steps)
- Epithelial cell states (healthy, injury levels 1-5)
- Phagocyte phenotypes (M0, M1, M2) and activation levels
- Microbe counts (commensals, pathogens)
- Treg states (resting, active)
- Cumulative death counts by mechanism
- Metadata: `param_set_id`, `scenario_id`, `replicate_id`, `seed`

## 🔧 Troubleshooting

### "Out of Memory" Error
- Reduce `N_PARAM_SETS` to 50 or fewer
- Reduce `N_REPLICATES` to 5
- Make sure `COMBINE_RESULTS = False`

### "Session Disconnected"
- This is normal for long-running simulations
- Re-run all cells in the notebook
- The simulation will skip already-completed parameter sets

### Slow Performance
- Make sure you've enabled compute units in your Colab settings
- Check the first cell output to see how many CPUs are available
- With compute units, you should have 8+ CPUs

### Import Errors
- Make sure you uploaded all 3 `*_G.py` files
- Re-run the package installation cell
- Restart the runtime if necessary

## 📈 Expected Performance

With compute units and proper setup:
- **Per simulation**: ~1-2 seconds
- **Per parameter set** (6 scenarios × 10 replicates): ~1-2 minutes
- **100 parameter sets**: ~2-3 hours

The actual time will vary based on your allocated resources.

## 💾 Downloading Results

### Option 1: Direct Download (Recommended)
The notebook creates a ZIP file of all results and downloads it to your computer.

### Option 2: Save to Google Drive
The notebook has an optional cell to copy results to your Google Drive for long-term storage.

## 🔄 Resuming Interrupted Simulations

If your simulation is interrupted:
1. Re-open the notebook
2. Re-upload the Python files (or mount Google Drive if you saved them there)
3. Use the **same** `OUTPUT_DIR` as before
4. Run the simulation cell again

The script will:
- Detect which parameter sets are already complete
- Skip them
- Continue with the remaining parameter sets

## 📚 Next Steps

After downloading your results:
1. Extract the ZIP file
2. Use your local environment (R or Python) for analysis
3. Combine parameter set files if needed:

```python
import pandas as pd
import glob

# Combine all parameter set files
files = glob.glob('mass_sim_results/simulation_results_param_set_*.csv')
all_results = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
all_results.to_csv('all_results_combined.csv', index=False)
```

## ❓ Questions?

If you encounter any issues:
1. Check the troubleshooting section above
2. Review the error message in the notebook
3. Make sure all dependencies are installed
4. Verify you have sufficient compute units

## 🆚 Differences from Local Version

| Feature | Local (`mass_simulation_LHS.py`) | Colab (`mass_simulation_LHS_G.py`) |
|---------|----------------------------------|-------------------------------------|
| Paths | Hardcoded absolute paths | Relative paths |
| Progress | Basic output | Detailed with ETA |
| CPU Detection | Manual setting | Auto-detect all available |
| Resume | Not supported | Automatic resume capability |
| File handling | System-specific | Cloud-optimized |

## 📝 License

These scripts are optimized versions of your original simulation code for Google Colab execution.
