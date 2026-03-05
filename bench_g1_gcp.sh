#!/usr/bin/env bash
# Benchmark G1 scalar multiplication on a GCP VM with AVX512-IFMA support.
# Prerequisites: gcloud CLI installed and authenticated (gcloud init / gcloud auth login).
#
# Usage: ./bench_g1_gcp.sh [--keep]
#   --keep   Don't delete the VM after benchmarking (for debugging)

set -euo pipefail

# ── Configuration ──────────────────────────────────────────────────────────
VM_NAME="vroom-bench"
ZONE="us-central1-a"
# C3 uses Intel Sapphire Rapids which has AVX512-IFMA
MACHINE_TYPE="c3-standard-4"
IMAGE_FAMILY="debian-12"
IMAGE_PROJECT="debian-cloud"
KEEP_VM=false

if [[ "${1:-}" == "--keep" ]]; then
    KEEP_VM=true
fi

# ── Preflight checks ──────────────────────────────────────────────────────
if ! command -v gcloud &>/dev/null; then
    echo "ERROR: gcloud CLI not found. Install it first:"
    echo "  brew install --cask google-cloud-sdk"
    exit 1
fi

PROJECT=$(gcloud config get-value project 2>/dev/null)
if [[ -z "$PROJECT" ]]; then
    echo "ERROR: No GCP project set. Run: gcloud init"
    exit 1
fi
echo "Using GCP project: $PROJECT"

# Check that Compute Engine API is enabled
if ! gcloud services list --enabled --filter="name:compute.googleapis.com" --format="value(name)" 2>/dev/null | grep -q compute; then
    echo "Enabling Compute Engine API..."
    gcloud services enable compute.googleapis.com
fi

# ── Create VM ──────────────────────────────────────────────────────────────
cleanup() {
    if [[ "$KEEP_VM" == false ]]; then
        echo ""
        echo "Cleaning up: deleting VM $VM_NAME..."
        gcloud compute instances delete "$VM_NAME" --zone="$ZONE" --quiet 2>/dev/null || true
    else
        echo ""
        echo "VM $VM_NAME kept alive. Delete manually when done:"
        echo "  gcloud compute instances delete $VM_NAME --zone=$ZONE --quiet"
    fi
}
trap cleanup EXIT

# Delete any leftover VM from a previous run and wait for the operation to complete
if gcloud compute instances describe "$VM_NAME" --zone="$ZONE" &>/dev/null; then
    echo "Deleting leftover VM from previous run..."
    gcloud compute instances delete "$VM_NAME" --zone="$ZONE" --quiet
    echo "Waiting for deletion to propagate..."
    sleep 15
fi

echo "Creating VM: $VM_NAME ($MACHINE_TYPE in $ZONE)..."
gcloud compute instances create "$VM_NAME" \
    --zone="$ZONE" \
    --machine-type="$MACHINE_TYPE" \
    --image-family="$IMAGE_FAMILY" \
    --image-project="$IMAGE_PROJECT" \
    --boot-disk-size=20GB \
    --boot-disk-type=pd-ssd \
    --quiet

echo "Waiting for VM to be ready..."
for i in $(seq 1 30); do
    if gcloud compute ssh "$VM_NAME" --zone="$ZONE" --command="echo ok" 2>/dev/null; then
        break
    fi
    if [[ $i -eq 30 ]]; then
        echo "ERROR: Timed out waiting for VM"
        exit 1
    fi
    sleep 5
done

# ── Install dependencies on VM ────────────────────────────────────────────
echo "Installing build dependencies..."
gcloud compute ssh "$VM_NAME" --zone="$ZONE" --command="
set -e
sudo apt-get update -qq
sudo apt-get install -y -qq clang llvm lld g++ libgmp-dev libbenchmark-dev make
"

# ── Upload source code ────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
echo "Uploading source code..."

# Create a tar excluding build artifacts and .git
tar czf /tmp/vroom-src.tar.gz \
    -C "$SCRIPT_DIR" \
    --exclude='.git' \
    --exclude='*.o' \
    --exclude='*.a' \
    --exclude='notes' \
    .

gcloud compute scp /tmp/vroom-src.tar.gz "$VM_NAME":~/vroom-src.tar.gz --zone="$ZONE" --quiet
rm /tmp/vroom-src.tar.gz

# ── Build and run benchmark ───────────────────────────────────────────────
echo "Building and running G1 scalar multiplication benchmark..."
gcloud compute ssh "$VM_NAME" --zone="$ZONE" --command="
set -e

# Verify AVX512-IFMA support
if grep -q avx512ifma /proc/cpuinfo; then
    echo 'AVX512-IFMA: supported'
else
    echo 'WARNING: AVX512-IFMA not detected on this machine!'
    cat /proc/cpuinfo | grep -m1 'model name'
fi

mkdir -p ~/vroom && cd ~/vroom
tar xzf ~/vroom-src.tar.gz

echo ''
echo '=== Building BLST ==='
cd ~/vroom/blst && make CC=clang

echo ''
echo '=== Building benchmark ==='
cd ~/vroom/src && make CXX=clang++ bench_pairing_50bit

echo ''
echo '=========================================='
echo '  G1 Scalar Multiplication Benchmark'
echo '=========================================='
~/vroom/src/bench_pairing_50bit --benchmark_filter='BM_G1_'
"

echo ""
echo "Benchmark complete."
