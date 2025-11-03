#!/bin/bash
################################################################################
# GCE Setup Script for Tregs Mass Simulation
#
# This script sets up a fresh Google Compute Engine instance to run
# mass_simulation_LHS.py with optimal configuration for 32-core machines.
#
# Usage:
#   chmod +x setup_gce.sh
#   ./setup_gce.sh
################################################################################

set -e  # Exit on any error

echo "========================================="
echo "GCE Setup for Tregs Simulation"
echo "========================================="
echo ""

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if running on Ubuntu
if [ ! -f /etc/lsb-release ]; then
    print_error "This script is designed for Ubuntu. Exiting."
    exit 1
fi

print_status "Detected OS: $(cat /etc/lsb-release | grep DESCRIPTION | cut -d'"' -f2)"

# Update system
print_status "Updating system packages..."
sudo apt-get update -qq
sudo apt-get upgrade -y -qq
print_success "System updated"

# Install essential tools
print_status "Installing essential tools (git, curl, wget, htop, tmux)..."
sudo apt-get install -y -qq git curl wget htop tmux build-essential
print_success "Essential tools installed"

# Install Python 3 and pip
print_status "Installing Python 3 and pip..."
sudo apt-get install -y -qq python3 python3-pip python3-venv python3-dev
PYTHON_VERSION=$(python3 --version)
print_success "Installed: $PYTHON_VERSION"

# Create project directory
print_status "Setting up project directory..."
cd ~
if [ ! -d "tregs" ]; then
    print_status "Cloning tregs repository..."
    read -p "Enter your GitHub username: " GITHUB_USER
    git clone https://github.com/$GITHUB_USER/tregs.git
    print_success "Repository cloned"
else
    print_warning "tregs directory already exists, skipping clone"
    cd tregs
    print_status "Pulling latest changes..."
    git pull
fi

cd ~/tregs

# Create virtual environment
print_status "Creating Python virtual environment..."
if [ ! -d "venv" ]; then
    python3 -m venv venv
    print_success "Virtual environment created"
else
    print_warning "Virtual environment already exists"
fi

# Activate virtual environment
print_status "Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
print_status "Upgrading pip..."
pip install --upgrade pip -q

# Install required Python packages
print_status "Installing Python packages (numpy, pandas, scipy, numba, matplotlib)..."
pip install numpy pandas scipy numba matplotlib -q
print_success "Python packages installed"

# Verify installations
print_status "Verifying installations..."

# Check Python packages
python3 << EOF
try:
    import numpy as np
    import pandas as pd
    import scipy
    import numba
    import matplotlib
    print(f"✓ numpy {np.__version__}")
    print(f"✓ pandas {pd.__version__}")
    print(f"✓ scipy {scipy.__version__}")
    print(f"✓ numba {numba.__version__}")
    print(f"✓ matplotlib {matplotlib.__version__}")
except ImportError as e:
    print(f"✗ Error importing package: {e}")
    exit(1)
EOF

# Check CPU cores
print_status "Checking system resources..."
NUM_CORES=$(python3 -c "import multiprocessing; print(multiprocessing.cpu_count())")
TOTAL_RAM=$(free -h | grep Mem | awk '{print $2}')

echo ""
echo "========================================="
echo "System Resources:"
echo "========================================="
echo "CPU Cores: $NUM_CORES"
echo "Total RAM: $TOTAL_RAM"
echo "Disk Space:"
df -h / | tail -1
echo "========================================="
echo ""

# Create helper scripts
print_status "Creating helper scripts..."

# Create run_simulation.sh
cat > run_simulation.sh << 'RUNSCRIPT'
#!/bin/bash
# Quick script to run simulations with common configurations

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPT_DIR

# Activate virtual environment
source venv/bin/activate

# Default values
N_PARAM_SETS=100
N_REPLICATES=10
N_CORES=$(python3 -c "import multiprocessing; print(multiprocessing.cpu_count())")
OUTPUT_DIR="mass_sim_results"
BASE_SEED=42

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --n_param_sets)
            N_PARAM_SETS="$2"
            shift 2
            ;;
        --n_replicates)
            N_REPLICATES="$2"
            shift 2
            ;;
        --n_cores)
            N_CORES="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --base_seed)
            BASE_SEED="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

echo "========================================"
echo "Running Mass Simulation"
echo "========================================"
echo "Parameter sets: $N_PARAM_SETS"
echo "Replicates: $N_REPLICATES"
echo "CPU cores: $N_CORES"
echo "Output dir: $OUTPUT_DIR"
echo "Base seed: $BASE_SEED"
echo "========================================"
echo ""

# Run simulation
python3 mass_simulation_LHS.py \
    --n_param_sets $N_PARAM_SETS \
    --n_replicates $N_REPLICATES \
    --n_cores $N_CORES \
    --output_dir $OUTPUT_DIR \
    --base_seed $BASE_SEED
RUNSCRIPT

chmod +x run_simulation.sh
print_success "Created run_simulation.sh"

# Create monitor_simulation.sh
cat > monitor_simulation.sh << 'MONITORSCRIPT'
#!/bin/bash
# Monitor simulation progress

OUTPUT_DIR="${1:-mass_sim_results}"

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Output directory $OUTPUT_DIR not found"
    exit 1
fi

echo "Monitoring simulation in: $OUTPUT_DIR"
echo "========================================"
echo ""

# Count completed parameter sets
COMPLETED=$(ls $OUTPUT_DIR/simulation_results_param_set_*.csv 2>/dev/null | wc -l)
echo "Completed parameter sets: $COMPLETED"

# Show disk usage
echo "Disk usage: $(du -sh $OUTPUT_DIR | cut -f1)"

# Show last 10 lines of newest file
NEWEST=$(ls -t $OUTPUT_DIR/simulation_results_param_set_*.csv 2>/dev/null | head -1)
if [ ! -z "$NEWEST" ]; then
    echo ""
    echo "Latest file: $(basename $NEWEST)"
    echo "Size: $(du -h $NEWEST | cut -f1)"
fi

echo ""
echo "========================================"
echo "System Resources:"
echo "========================================"
echo "CPU Usage:"
top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1"%"}'

echo "Memory Usage:"
free -h | grep Mem | awk '{print $3 " / " $2}'

echo "Disk Usage:"
df -h / | tail -1 | awk '{print $3 " / " $2 " (" $5 ")"}'
echo "========================================"
MONITORSCRIPT

chmod +x monitor_simulation.sh
print_success "Created monitor_simulation.sh"

# Create test_setup.sh
cat > test_setup.sh << 'TESTSCRIPT'
#!/bin/bash
# Quick test to verify everything works

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $SCRIPT_DIR

# Activate virtual environment
source venv/bin/activate

echo "Running quick test simulation..."
echo "This will run 2 parameter sets with 2 replicates"
echo ""

python3 mass_simulation_LHS.py \
    --n_param_sets 2 \
    --n_replicates 2 \
    --n_cores 4 \
    --output_dir test_results

if [ $? -eq 0 ]; then
    echo ""
    echo "✓ Test successful!"
    echo "Check test_results/ for output files"
else
    echo ""
    echo "✗ Test failed"
    exit 1
fi
TESTSCRIPT

chmod +x test_setup.sh
print_success "Created test_setup.sh"

# Add helpful aliases to .bashrc
print_status "Adding helpful aliases to ~/.bashrc..."
cat >> ~/.bashrc << 'ALIASES'

# Tregs simulation aliases
alias tregs='cd ~/tregs && source venv/bin/activate'
alias run_sim='cd ~/tregs && source venv/bin/activate && ./run_simulation.sh'
alias monitor_sim='cd ~/tregs && ./monitor_simulation.sh'
alias test_sim='cd ~/tregs && ./test_setup.sh'
ALIASES

print_success "Aliases added (run 'source ~/.bashrc' to activate)"

# Create README for quick reference
cat > ~/tregs/GCE_QUICKSTART.md << 'QUICKSTART'
# GCE Quick Start

Your instance is now set up! Here's how to use it:

## Quick Commands

```bash
# Activate environment (or just type 'tregs')
cd ~/tregs
source venv/bin/activate

# Test the setup (2 minutes)
./test_setup.sh

# Run full simulation (use all CPU cores)
./run_simulation.sh --n_param_sets 100 --n_replicates 10

# Monitor progress
./monitor_simulation.sh

# Run in background with tmux
tmux new -s sim
./run_simulation.sh --n_param_sets 1000 --n_replicates 100
# Press Ctrl+B, then D to detach
# Reattach: tmux attach -t sim
```

## Custom Configuration

```bash
./run_simulation.sh \
    --n_param_sets 500 \
    --n_replicates 50 \
    --n_cores 32 \
    --output_dir my_results \
    --base_seed 12345
```

## Download Results (from your local computer)

```bash
# Download entire directory
gcloud compute scp --recurse \
    INSTANCE_NAME:~/tregs/mass_sim_results \
    ./local_results \
    --zone=YOUR_ZONE

# Compress first (recommended for large results)
# On GCE:
tar -czf results.tar.gz mass_sim_results/

# On local:
gcloud compute scp INSTANCE_NAME:~/tregs/results.tar.gz ./ --zone=YOUR_ZONE
```

## System Monitoring

```bash
# Real-time CPU/memory monitor
htop

# Check disk space
df -h

# Count completed parameter sets
ls mass_sim_results/simulation_results_param_set_*.csv | wc -l
```

## Helpful Aliases

After running `source ~/.bashrc`, you can use:
- `tregs` - cd to project and activate venv
- `run_sim` - quick run with default params
- `monitor_sim` - check progress
- `test_sim` - run quick test

## Next Steps

1. Run test: `./test_setup.sh`
2. Run full simulation: `./run_simulation.sh`
3. Download results when done
4. Stop instance to save costs: `gcloud compute instances stop INSTANCE_NAME --zone=YOUR_ZONE`

See GCE_SETUP_GUIDE.md for full documentation.
QUICKSTART

print_success "Created GCE_QUICKSTART.md"

# Final summary
echo ""
echo "========================================="
echo "Setup Complete!"
echo "========================================="
echo ""
print_success "Your GCE instance is ready to run simulations!"
echo ""
echo "Next steps:"
echo "  1. Test the setup: ./test_setup.sh"
echo "  2. Run simulation: ./run_simulation.sh --n_param_sets 100 --n_replicates 10"
echo "  3. Monitor progress: ./monitor_simulation.sh"
echo ""
echo "Quick reference: cat GCE_QUICKSTART.md"
echo "Full guide: cat GCE_SETUP_GUIDE.md"
echo ""
echo "Activate environment: source venv/bin/activate"
echo "Or reload aliases: source ~/.bashrc"
echo "Then just type: tregs"
echo ""
echo "========================================="
echo "System is ready with $NUM_CORES CPU cores and $TOTAL_RAM RAM"
echo "========================================="
