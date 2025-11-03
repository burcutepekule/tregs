# GCE Quick Start Guide

This is a quick reference for setting up and running simulations on Google Compute Engine with 32 cores.

## Option 1: Automated Setup (Recommended)

### Step 1: Create Instance (from your local computer)

```bash
# Download the creation script
curl -O https://raw.githubusercontent.com/burcutepekule/tregs/main/create_gce_instance.sh

# Make it executable
chmod +x create_gce_instance.sh

# Run it (follow the prompts)
./create_gce_instance.sh
```

This will:
- Create a 32-core GCE instance
- Open SSH connection automatically
- Save connection details to `gce_instance_info.txt`

### Step 2: Setup the Instance (once SSH connected)

```bash
# Clone repository
git clone https://github.com/burcutepekule/tregs.git
cd tregs

# Run automated setup
chmod +x setup_gce.sh
./setup_gce.sh
```

This will:
- Install Python and dependencies
- Set up virtual environment
- Install numpy, pandas, scipy, numba
- Create helper scripts
- Verify everything works

### Step 3: Run Simulation

```bash
# Quick test (2 minutes)
./test_setup.sh

# Full production run (100 parameter sets, 10 replicates)
./run_simulation.sh --n_param_sets 100 --n_replicates 10

# Large run (1000 parameter sets, 100 replicates, 8-12 hours)
./run_simulation.sh --n_param_sets 1000 --n_replicates 100
```

### Step 4: Download Results (from your local computer)

```bash
# Compress on GCE
gcloud compute ssh YOUR_INSTANCE_NAME --zone=YOUR_ZONE --command="cd ~/tregs && tar -czf results.tar.gz mass_sim_results/"

# Download
gcloud compute scp YOUR_INSTANCE_NAME:~/tregs/results.tar.gz ./ --zone=YOUR_ZONE

# Extract
tar -xzf results.tar.gz
```

---

## Option 2: Manual gcloud Commands

### Create Instance

```bash
gcloud compute instances create tregs-simulation-32core \
    --zone=us-central1-a \
    --machine-type=n2-standard-32 \
    --boot-disk-size=50GB \
    --boot-disk-type=pd-balanced \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud
```

### Connect

```bash
gcloud compute ssh tregs-simulation-32core --zone=us-central1-a
```

### Setup (on the instance)

```bash
# Update system
sudo apt-get update && sudo apt-get upgrade -y

# Install dependencies
sudo apt-get install -y python3 python3-pip python3-venv git htop tmux

# Clone repo
git clone https://github.com/burcutepekule/tregs.git
cd tregs

# Setup Python environment
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install numpy pandas scipy numba matplotlib
```

### Run Simulation

```bash
python3 mass_simulation_LHS.py \
    --n_param_sets 100 \
    --n_replicates 10 \
    --n_cores 32 \
    --output_dir mass_sim_results
```

---

## Using tmux for Long-Running Jobs

```bash
# Start tmux session
tmux new -s simulation

# Run your simulation
./run_simulation.sh --n_param_sets 1000 --n_replicates 100

# Detach from session: Press Ctrl+B, then D
# You can now close your SSH connection safely

# Reattach later
gcloud compute ssh YOUR_INSTANCE --zone=YOUR_ZONE
tmux attach -t simulation

# List sessions
tmux ls

# Kill session (after job completes)
tmux kill-session -t simulation
```

---

## Monitoring

### Check Progress

```bash
# Use the helper script
./monitor_simulation.sh

# Or manually:
# Count completed parameter sets
ls mass_sim_results/simulation_results_param_set_*.csv | wc -l

# Check disk usage
du -sh mass_sim_results/

# Monitor CPU/memory
htop
```

### Check if Simulation is Running

```bash
# Check Python processes
ps aux | grep mass_simulation

# Check from tmux
tmux ls
tmux attach -t simulation
```

---

## Cost Management

### Estimated Costs (n2-standard-32: 32 vCPU, 128GB RAM)

- **On-demand**: ~$1.54/hour (~$37/day)
- **Spot VM**: ~$0.46/hour (~$11/day)

### Stop Instance When Not in Use

```bash
# From local computer
gcloud compute instances stop tregs-simulation-32core --zone=us-central1-a

# Start again
gcloud compute instances start tregs-simulation-32core --zone=us-central1-a
```

**Note**: Stopped instances still incur storage costs (~$2/month for 50GB), but save 90%+ on compute costs.

### Delete Instance (When Completely Done)

```bash
# CAUTION: This deletes everything!
gcloud compute instances delete tregs-simulation-32core --zone=us-central1-a
```

---

## Example Workflows

### Workflow 1: Quick Test

```bash
# Create instance
./create_gce_instance.sh

# On instance: setup and test
./setup_gce.sh
./test_setup.sh

# If successful, stop instance
exit
gcloud compute instances stop YOUR_INSTANCE --zone=YOUR_ZONE
```

### Workflow 2: Production Run

```bash
# Start instance (if stopped)
gcloud compute instances start YOUR_INSTANCE --zone=YOUR_ZONE

# SSH in
gcloud compute ssh YOUR_INSTANCE --zone=YOUR_ZONE

# Run in tmux
tmux new -s sim
source venv/bin/activate
./run_simulation.sh --n_param_sets 1000 --n_replicates 100

# Detach: Ctrl+B, then D
# Close SSH (job keeps running)

# Check later
gcloud compute ssh YOUR_INSTANCE --zone=YOUR_ZONE
tmux attach -t sim

# When done, download results
exit
gcloud compute scp --recurse YOUR_INSTANCE:~/tregs/mass_sim_results ./ --zone=YOUR_ZONE

# Stop instance
gcloud compute instances stop YOUR_INSTANCE --zone=YOUR_ZONE
```

### Workflow 3: Multiple Runs in Parallel

```bash
# SSH into instance
cd ~/tregs
source venv/bin/activate

# Terminal 1: First batch
tmux new -s batch1
./run_simulation.sh --n_param_sets 500 --output_dir batch1 --base_seed 1
# Ctrl+B, D

# Terminal 2: Second batch
tmux new -s batch2
./run_simulation.sh --n_param_sets 500 --output_dir batch2 --base_seed 1000
# Ctrl+B, D

# Monitor both
tmux attach -t batch1  # Check, then Ctrl+B, D
tmux attach -t batch2  # Check, then Ctrl+B, D
```

---

## Troubleshooting

### Issue: Can't SSH into instance

```bash
# Check if instance is running
gcloud compute instances list

# Start if stopped
gcloud compute instances start YOUR_INSTANCE --zone=YOUR_ZONE

# Try SSH via browser
# Go to console.cloud.google.com → Compute Engine → VM Instances → Click SSH
```

### Issue: Out of memory

```bash
# Check memory
free -h

# Reduce simulation size
./run_simulation.sh --n_param_sets 50 --n_replicates 5
```

### Issue: Simulation seems stuck

```bash
# Check if running
ps aux | grep python

# Check CPU usage (should be ~3200% for 32 cores)
htop

# Check logs
tail -f simulation.log  # if using nohup
```

### Issue: Not using all CPU cores

```bash
# Verify core count
python3 -c "import multiprocessing; print(multiprocessing.cpu_count())"

# Make sure you specified --n_cores 32
./run_simulation.sh --n_cores 32
```

---

## Quick Reference Card

```bash
# CREATE INSTANCE (local)
./create_gce_instance.sh

# SETUP (on GCE)
git clone https://github.com/burcutepekule/tregs.git
cd tregs
./setup_gce.sh

# TEST (on GCE)
./test_setup.sh

# RUN (on GCE)
./run_simulation.sh --n_param_sets 100 --n_replicates 10

# MONITOR (on GCE)
./monitor_simulation.sh

# TMUX (on GCE)
tmux new -s sim          # Start
./run_simulation.sh      # Run job
# Ctrl+B, D              # Detach
tmux attach -t sim       # Reattach

# DOWNLOAD (local)
gcloud compute scp --recurse \
  YOUR_INSTANCE:~/tregs/mass_sim_results ./ \
  --zone=YOUR_ZONE

# STOP (local)
gcloud compute instances stop YOUR_INSTANCE --zone=YOUR_ZONE

# START (local)
gcloud compute instances start YOUR_INSTANCE --zone=YOUR_ZONE

# DELETE (local)
gcloud compute instances delete YOUR_INSTANCE --zone=YOUR_ZONE
```

---

## Performance Estimates

With 32 cores:

| Simulation Size | Time |
|----------------|------|
| Test (2 param sets, 2 reps) | 1-2 min |
| Small (10 param sets, 10 reps) | 5-10 min |
| Medium (100 param sets, 10 reps) | 30-60 min |
| Large (100 param sets, 100 reps) | 5-8 hours |
| Very Large (1000 param sets, 100 reps) | 2-3 days |

**Note**: Times are approximate and depend on specific parameters and system load.

---

## Getting Help

- Full documentation: `cat GCE_SETUP_GUIDE.md`
- GCP documentation: https://cloud.google.com/compute/docs
- Check instance details: `cat gce_instance_info.txt`

---

## After You're Done

1. Download all results
2. Stop or delete the instance to avoid charges
3. Verify results are safe locally
4. Optional: Upload to Cloud Storage for backup

```bash
# Backup to Cloud Storage (optional)
gsutil mb gs://your-bucket-name-tregs-results
gsutil -m cp -r mass_sim_results/ gs://your-bucket-name-tregs-results/
```
