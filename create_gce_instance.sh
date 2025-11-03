#!/bin/bash
################################################################################
# Create GCE Instance Script for Tregs Simulation
#
# This script creates a GCE instance with 32 cores optimized for running
# mass_simulation_LHS.py simulations.
#
# Prerequisites:
#   - gcloud CLI installed and configured
#   - GCP project with billing enabled
#   - Compute Engine API enabled
#
# Usage:
#   chmod +x create_gce_instance.sh
#   ./create_gce_instance.sh
################################################################################

set -e

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

print_status() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

echo "========================================="
echo "GCE Instance Creator for Tregs Simulation"
echo "========================================="
echo ""

# Check if gcloud is installed
if ! command -v gcloud &> /dev/null; then
    print_error "gcloud CLI is not installed"
    echo "Please install it from: https://cloud.google.com/sdk/docs/install"
    exit 1
fi

print_success "gcloud CLI found"

# Get current project
CURRENT_PROJECT=$(gcloud config get-value project 2>/dev/null)
if [ -z "$CURRENT_PROJECT" ]; then
    print_error "No GCP project is configured"
    echo "Run: gcloud config set project YOUR_PROJECT_ID"
    exit 1
fi

print_status "Current GCP project: $CURRENT_PROJECT"

# Configuration
echo ""
echo "========================================="
echo "Instance Configuration"
echo "========================================="
echo ""

# Instance name
read -p "Instance name [tregs-simulation-32core]: " INSTANCE_NAME
INSTANCE_NAME=${INSTANCE_NAME:-tregs-simulation-32core}

# Zone
echo ""
echo "Recommended zones:"
echo "  us-central1-a (Iowa)"
echo "  us-east1-b (South Carolina)"
echo "  us-west1-b (Oregon)"
echo "  europe-west1-b (Belgium)"
echo "  asia-east1-b (Taiwan)"
read -p "Zone [us-central1-a]: " ZONE
ZONE=${ZONE:-us-central1-a}

# Machine type
echo ""
echo "Available machine types:"
echo "  1. n2-standard-32   (32 vCPU, 128 GB RAM) - Balanced [RECOMMENDED] ~\$1.54/hr"
echo "  2. n2-highcpu-32    (32 vCPU, 32 GB RAM)  - CPU optimized ~\$1.16/hr"
echo "  3. n2-highmem-32    (32 vCPU, 256 GB RAM) - Memory optimized ~\$2.32/hr"
echo "  4. c2-standard-60   (60 vCPU, 240 GB RAM) - Compute optimized ~\$3.10/hr"
echo "  5. n2-standard-64   (64 vCPU, 256 GB RAM) - More cores ~\$3.08/hr"
echo "  6. Custom"
read -p "Choose machine type [1]: " MACHINE_CHOICE
MACHINE_CHOICE=${MACHINE_CHOICE:-1}

case $MACHINE_CHOICE in
    1) MACHINE_TYPE="n2-standard-32" ;;
    2) MACHINE_TYPE="n2-highcpu-32" ;;
    3) MACHINE_TYPE="n2-highmem-32" ;;
    4) MACHINE_TYPE="c2-standard-60" ;;
    5) MACHINE_TYPE="n2-standard-64" ;;
    6)
        read -p "Enter custom machine type: " MACHINE_TYPE
        ;;
    *)
        print_warning "Invalid choice, using n2-standard-32"
        MACHINE_TYPE="n2-standard-32"
        ;;
esac

# Preemptible/Spot option
echo ""
read -p "Use Spot/Preemptible VM? (70% cheaper, can be interrupted) [y/N]: " USE_PREEMPTIBLE
PREEMPTIBLE_FLAG=""
if [[ $USE_PREEMPTIBLE =~ ^[Yy]$ ]]; then
    PREEMPTIBLE_FLAG="--provisioning-model=SPOT --instance-termination-action=STOP"
    print_warning "Using Spot VM - instance may be interrupted!"
fi

# Boot disk size
read -p "Boot disk size in GB [50]: " DISK_SIZE
DISK_SIZE=${DISK_SIZE:-50}

# Confirm configuration
echo ""
echo "========================================="
echo "Configuration Summary"
echo "========================================="
echo "Instance name: $INSTANCE_NAME"
echo "Zone: $ZONE"
echo "Machine type: $MACHINE_TYPE"
echo "Spot/Preemptible: ${USE_PREEMPTIBLE:-No}"
echo "Disk size: ${DISK_SIZE}GB"
echo "Project: $CURRENT_PROJECT"
echo "========================================="
echo ""
read -p "Create instance with this configuration? [Y/n]: " CONFIRM
if [[ $CONFIRM =~ ^[Nn]$ ]]; then
    print_warning "Cancelled"
    exit 0
fi

# Create instance
print_status "Creating instance '$INSTANCE_NAME'..."
echo ""

gcloud compute instances create $INSTANCE_NAME \
    --zone=$ZONE \
    --machine-type=$MACHINE_TYPE \
    --boot-disk-size=${DISK_SIZE}GB \
    --boot-disk-type=pd-balanced \
    --image-family=ubuntu-2204-lts \
    --image-project=ubuntu-os-cloud \
    --metadata=startup-script='#!/bin/bash
apt-get update
apt-get install -y git python3-pip python3-venv' \
    $PREEMPTIBLE_FLAG

if [ $? -eq 0 ]; then
    echo ""
    print_success "Instance created successfully!"
    echo ""
    echo "========================================="
    echo "Next Steps"
    echo "========================================="
    echo ""
    echo "1. SSH into your instance:"
    echo "   gcloud compute ssh $INSTANCE_NAME --zone=$ZONE"
    echo ""
    echo "2. Download and run setup script:"
    echo "   wget https://raw.githubusercontent.com/YOUR_USERNAME/tregs/main/setup_gce.sh"
    echo "   chmod +x setup_gce.sh"
    echo "   ./setup_gce.sh"
    echo ""
    echo "   OR manually:"
    echo "   git clone https://github.com/YOUR_USERNAME/tregs.git"
    echo "   cd tregs"
    echo "   chmod +x setup_gce.sh"
    echo "   ./setup_gce.sh"
    echo ""
    echo "3. After setup, run test:"
    echo "   ./test_setup.sh"
    echo ""
    echo "4. Run your simulation:"
    echo "   ./run_simulation.sh --n_param_sets 100 --n_replicates 10"
    echo ""
    echo "========================================="
    echo "Management Commands"
    echo "========================================="
    echo ""
    echo "SSH into instance:"
    echo "  gcloud compute ssh $INSTANCE_NAME --zone=$ZONE"
    echo ""
    echo "Stop instance (to save costs):"
    echo "  gcloud compute instances stop $INSTANCE_NAME --zone=$ZONE"
    echo ""
    echo "Start instance:"
    echo "  gcloud compute instances start $INSTANCE_NAME --zone=$ZONE"
    echo ""
    echo "Delete instance (CAUTION - deletes everything!):"
    echo "  gcloud compute instances delete $INSTANCE_NAME --zone=$ZONE"
    echo ""
    echo "Download results (from local computer):"
    echo "  gcloud compute scp --recurse $INSTANCE_NAME:~/tregs/mass_sim_results ./ --zone=$ZONE"
    echo ""
    echo "========================================="
    echo ""

    # Create a local reference file
    cat > gce_instance_info.txt << EOF
GCE Instance Information
========================

Instance Name: $INSTANCE_NAME
Zone: $ZONE
Machine Type: $MACHINE_TYPE
Spot/Preemptible: ${USE_PREEMPTIBLE:-No}
Project: $CURRENT_PROJECT
Created: $(date)

Quick Commands:
---------------

SSH:
  gcloud compute ssh $INSTANCE_NAME --zone=$ZONE

Stop:
  gcloud compute instances stop $INSTANCE_NAME --zone=$ZONE

Start:
  gcloud compute instances start $INSTANCE_NAME --zone=$ZONE

Delete:
  gcloud compute instances delete $INSTANCE_NAME --zone=$ZONE

Download Results:
  gcloud compute scp --recurse $INSTANCE_NAME:~/tregs/mass_sim_results ./ --zone=$ZONE

View in Console:
  https://console.cloud.google.com/compute/instances?project=$CURRENT_PROJECT
EOF

    print_success "Instance info saved to: gce_instance_info.txt"
    echo ""
    print_status "Opening SSH connection in 5 seconds..."
    sleep 5
    gcloud compute ssh $INSTANCE_NAME --zone=$ZONE
else
    print_error "Failed to create instance"
    echo ""
    echo "Common issues:"
    echo "  - Compute Engine API not enabled"
    echo "  - Insufficient quota for the machine type"
    echo "  - Billing not enabled on project"
    echo ""
    echo "Check the error message above for details."
    exit 1
fi
