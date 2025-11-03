# Running Mass Simulations on Google Cloud Compute Engine (32 Cores)

This guide shows you how to set up a Google Cloud Compute Engine (GCE) instance with 32 cores to run your `mass_simulation_LHS.py` simulations at maximum performance.

## Why GCE instead of Colab?

| Feature | Google Colab | GCE (32 cores) |
|---------|--------------|----------------|
| CPUs | 2-8 cores | 32+ cores (your choice) |
| RAM | 12-25 GB | Up to 208 GB |
| Session timeout | 12 hours max | No timeout (runs continuously) |
| Cost | Free/subscription | Pay per use (~$1.20/hour for 32 cores) |
| Stability | Can disconnect | Stable SSH connection |
| Performance | Limited | 4-10x faster with 32 cores |

## Prerequisites

1. A Google Cloud Platform (GCP) account
2. Billing enabled on your GCP project
3. Basic familiarity with terminal/SSH

## Part 1: Create GCE Instance with 32 Cores

### Step 1: Set Up Google Cloud Project

1. Go to [Google Cloud Console](https://console.cloud.google.com/)
2. Create a new project or select an existing one
3. Enable the Compute Engine API:
   - Navigate to "APIs & Services" → "Library"
   - Search for "Compute Engine API"
   - Click "Enable"

### Step 2: Create VM Instance

**Option A: Using the Web Console (Recommended for beginners)**

1. Go to **Compute Engine** → **VM Instances**
2. Click **"CREATE INSTANCE"**
3. Configure your instance:

   **Basic Configuration:**
   - **Name**: `tregs-simulation-32core` (or your choice)
   - **Region**: Choose one close to you (e.g., `us-central1`)
   - **Zone**: Any zone (e.g., `us-central1-a`)

   **Machine Configuration:**
   - **Machine family**: General-purpose
   - **Series**: N2 or N2D (best performance for compute)
   - **Machine type**: Click "CUSTOM" or select:
     - **n2-standard-32** (32 vCPUs, 128 GB RAM) - Recommended
     - **n2-highcpu-32** (32 vCPUs, 32 GB RAM) - Lower cost
     - **n2-highmem-32** (32 vCPUs, 256 GB RAM) - If you need more memory

   **Boot Disk:**
   - Click "CHANGE"
   - **Operating System**: Ubuntu
   - **Version**: Ubuntu 22.04 LTS
   - **Boot disk type**: Balanced persistent disk
   - **Size**: 50 GB (sufficient for your simulations)

   **Firewall:**
   - ✓ Allow HTTP traffic (optional)
   - ✓ Allow HTTPS traffic (optional)

4. Click **"CREATE"**

**Option B: Using gcloud CLI (Advanced)**

```bash
gcloud compute instances create tregs-simulation-32core \
    --zone=us-central1-a \
    --machine-type=n2-standard-32 \
    --boot-disk-size=50GB \
    --boot-disk-type=pd-balanced \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud
```

### Step 3: Connect to Your Instance

**Via Web Console:**
1. Go to **VM Instances**
2. Click **SSH** next to your instance name

**Via Terminal (Mac/Linux):**
```bash
gcloud compute ssh tregs-simulation-32core --zone=us-central1-a
```

**Via Terminal (Windows):**
Use the Cloud Shell or install gcloud CLI, then use the same command as above.

## Part 2: Set Up the Instance (Automated)

Once connected via SSH, run the automated setup script:

```bash
# Download the setup script
curl -O https://raw.githubusercontent.com/YOUR_USERNAME/tregs/main/setup_gce.sh

# Make it executable
chmod +x setup_gce.sh

# Run it
./setup_gce.sh
```

This script will:
1. Update Ubuntu packages
2. Install Python 3.10+
3. Install required Python packages (numpy, pandas, scipy, numba)
4. Clone your repository
5. Verify the setup

**OR manually follow Part 3 below if you prefer step-by-step setup.**

## Part 3: Manual Setup (Alternative to Automated)

### Step 1: Update System

```bash
sudo apt-get update
sudo apt-get upgrade -y
```

### Step 2: Install Python and Dependencies

```bash
# Install Python 3.10+
sudo apt-get install -y python3 python3-pip python3-venv

# Verify Python version
python3 --version  # Should be 3.10 or higher
```

### Step 3: Clone Your Repository

```bash
# Install git if not present
sudo apt-get install -y git

# Clone your repository
git clone https://github.com/YOUR_USERNAME/tregs.git
cd tregs
```

### Step 4: Set Up Python Environment

```bash
# Create virtual environment
python3 -m venv venv

# Activate it
source venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install required packages
pip install numpy pandas scipy numba matplotlib
```

### Step 5: Verify Setup

```bash
# Check CPU cores
python3 -c "import multiprocessing; print(f'Available CPUs: {multiprocessing.cpu_count()}')"
# Should show: Available CPUs: 32

# Test import
python3 -c "import numpy, pandas, scipy, numba; print('All packages imported successfully!')"
```

## Part 4: Run Your Simulation

### Option 1: Quick Test (2 minutes)

```bash
cd ~/tregs
source venv/bin/activate  # If not already activated

# Run a small test
python3 mass_simulation_LHS.py \
    --n_param_sets 2 \
    --n_replicates 2 \
    --n_cores 32 \
    --output_dir test_results
```

### Option 2: Full Production Run

```bash
# Example: 100 parameter sets, 10 replicates per scenario
python3 mass_simulation_LHS.py \
    --n_param_sets 100 \
    --n_replicates 10 \
    --n_cores 32 \
    --output_dir mass_sim_results \
    --base_seed 42

# Total simulations: 100 param sets × 6 scenarios × 10 replicates = 6,000 simulations
# Estimated time with 32 cores: ~30-60 minutes
```

### Option 3: Large-Scale Run (Overnight)

```bash
# Example: 1000 parameter sets, 100 replicates per scenario
python3 mass_simulation_LHS.py \
    --n_param_sets 1000 \
    --n_replicates 100 \
    --n_cores 32 \
    --output_dir large_sim_results \
    --base_seed 42

# Total simulations: 1000 × 6 × 100 = 600,000 simulations
# Estimated time with 32 cores: ~8-12 hours
```

### Run in Background (Recommended for Long Jobs)

Use `tmux` or `screen` to keep the job running even if you disconnect:

```bash
# Install tmux
sudo apt-get install -y tmux

# Start a tmux session
tmux new -s simulation

# Run your simulation
python3 mass_simulation_LHS.py --n_param_sets 1000 --n_replicates 100 --n_cores 32

# Detach from session: Press Ctrl+B, then D
# Reattach later: tmux attach -t simulation
# Check if running: tmux ls
```

Or use `nohup` to run in background:

```bash
nohup python3 mass_simulation_LHS.py \
    --n_param_sets 1000 \
    --n_replicates 100 \
    --n_cores 32 \
    --output_dir mass_sim_results \
    > simulation.log 2>&1 &

# Check progress
tail -f simulation.log

# Check if still running
ps aux | grep mass_simulation
```

## Part 5: Monitor Your Simulation

### Check CPU Usage

```bash
# Install htop for better monitoring
sudo apt-get install -y htop

# Monitor in real-time
htop
# Press F10 to quit
```

### Check Disk Usage

```bash
# Check disk space
df -h

# Check size of output directory
du -sh mass_sim_results/
```

### Check Progress

```bash
# Count completed parameter sets
ls mass_sim_results/simulation_results_param_set_*.csv | wc -l

# View the log
tail -f simulation.log  # If using nohup
```

## Part 6: Download Results

### Option 1: Use gcloud CLI (Recommended)

On your **local computer**, run:

```bash
# Download entire results directory
gcloud compute scp --recurse \
    tregs-simulation-32core:~/tregs/mass_sim_results \
    ./local_results \
    --zone=us-central1-a

# Download specific files only
gcloud compute scp \
    tregs-simulation-32core:~/tregs/mass_sim_results/*.csv \
    ./local_results/ \
    --zone=us-central1-a
```

### Option 2: Compress and Download

On the **GCE instance**:

```bash
# Compress results
cd ~/tregs
tar -czf mass_sim_results.tar.gz mass_sim_results/

# Download from local computer
gcloud compute scp \
    tregs-simulation-32core:~/tregs/mass_sim_results.tar.gz \
    ./ \
    --zone=us-central1-a

# Extract on local computer
tar -xzf mass_sim_results.tar.gz
```

### Option 3: Upload to Google Cloud Storage

```bash
# Install gsutil (usually pre-installed)
which gsutil

# Create a bucket (one-time setup)
gsutil mb gs://your-bucket-name-tregs-results

# Upload results
gsutil -m cp -r mass_sim_results/ gs://your-bucket-name-tregs-results/

# Download from anywhere
gsutil -m cp -r gs://your-bucket-name-tregs-results/mass_sim_results/ ./
```

## Part 7: Cost Management

### Estimate Costs

**n2-standard-32 (32 vCPUs, 128 GB RAM):**
- **On-demand**: ~$1.54/hour (~$37/day)
- **Preemptible**: ~$0.46/hour (~$11/day) - Can be interrupted!

**n2-highcpu-32 (32 vCPUs, 32 GB RAM):**
- **On-demand**: ~$1.16/hour (~$28/day)
- **Preemptible**: ~$0.35/hour (~$8/day)

### Save Money

1. **Use Preemptible/Spot VMs** (70% cheaper):
   - Good for jobs that can resume/restart
   - Your script already supports resume (skips completed param sets)
   - Add `--preemptible` flag when creating instance

2. **Stop instance when not in use**:
   ```bash
   # Stop (you'll still pay for disk, but not CPU)
   gcloud compute instances stop tregs-simulation-32core --zone=us-central1-a

   # Start again
   gcloud compute instances start tregs-simulation-32core --zone=us-central1-a
   ```

3. **Delete instance when completely done**:
   ```bash
   # CAUTION: This deletes everything!
   gcloud compute instances delete tregs-simulation-32core --zone=us-central1-a
   ```

4. **Set up budget alerts**:
   - Go to "Billing" → "Budgets & alerts"
   - Set alerts at $50, $100, etc.

### Monitor Costs

- Check real-time costs: [Billing Dashboard](https://console.cloud.google.com/billing)
- View detailed breakdown: "Billing" → "Reports"

## Part 8: Optimization Tips

### Maximize Performance

1. **Use all 32 cores**:
   ```bash
   python3 mass_simulation_LHS.py --n_cores 32
   ```

2. **Run multiple jobs in parallel** (if you have multiple scenarios):
   ```bash
   # Terminal 1
   python3 mass_simulation_LHS.py --n_param_sets 500 --output_dir batch1 --base_seed 1 &

   # Terminal 2
   python3 mass_simulation_LHS.py --n_param_sets 500 --output_dir batch2 --base_seed 1000 &
   ```

3. **Monitor memory usage**:
   ```bash
   free -h
   # If memory is low, reduce n_replicates or n_param_sets
   ```

### Resume Interrupted Jobs

If your simulation is interrupted:

1. **Restart the instance** (if stopped)
2. **SSH back in**
3. **Run the same command** - it will skip completed parameter sets:
   ```bash
   python3 mass_simulation_LHS.py --n_param_sets 1000 --n_cores 32 --output_dir mass_sim_results
   # Automatically resumes!
   ```

## Part 9: Troubleshooting

### Problem: "Permission denied" errors

```bash
# Make sure you're in your home directory
cd ~

# Check file permissions
ls -la tregs/
```

### Problem: "Out of memory"

```bash
# Check available memory
free -h

# Reduce batch size
python3 mass_simulation_LHS.py --n_param_sets 50 --n_replicates 5
```

### Problem: Instance is slow

```bash
# Verify you got 32 cores
nproc
# Should show: 32

# Check if CPU is being used
htop
```

### Problem: Can't connect via SSH

1. Check if instance is running in GCP Console
2. Check firewall rules allow SSH (port 22)
3. Try browser-based SSH from Console

### Problem: Simulation crashes

```bash
# Check error logs
cat simulation.log

# Test with small parameters first
python3 mass_simulation_LHS.py --n_param_sets 1 --n_replicates 1 --n_cores 2
```

## Part 10: Next Steps

### After Simulation Completes

1. **Download results** (see Part 6)
2. **Analyze locally**:
   ```bash
   python3 analyze_mass_sim.py --results_dir mass_sim_results
   ```
3. **Stop or delete instance** to save costs
4. **Backup to Cloud Storage** for long-term retention

### Scaling Up

If you need even more power:

- **n2-standard-64** (64 vCPUs, 256 GB RAM)
- **n2-standard-96** (96 vCPUs, 384 GB RAM)
- Multiple instances running in parallel

### Alternative: Using Compute-Optimized VMs

For pure CPU work, consider:
- **c2-standard-60** (60 vCPUs, 240 GB RAM)
- Higher per-core performance than N2

## Quick Reference Commands

```bash
# Create instance
gcloud compute instances create tregs-simulation-32core \
    --zone=us-central1-a --machine-type=n2-standard-32

# SSH into instance
gcloud compute ssh tregs-simulation-32core --zone=us-central1-a

# Run simulation in background
tmux new -s sim
python3 mass_simulation_LHS.py --n_param_sets 100 --n_cores 32
# Ctrl+B, then D to detach

# Reattach to session
tmux attach -t sim

# Download results
gcloud compute scp --recurse \
    tregs-simulation-32core:~/tregs/mass_sim_results ./

# Stop instance
gcloud compute instances stop tregs-simulation-32core --zone=us-central1-a

# Delete instance (CAUTION!)
gcloud compute instances delete tregs-simulation-32core --zone=us-central1-a
```

## Support

- [GCE Documentation](https://cloud.google.com/compute/docs)
- [gcloud CLI Reference](https://cloud.google.com/sdk/gcloud/reference/compute)
- [GCP Pricing Calculator](https://cloud.google.com/products/calculator)
