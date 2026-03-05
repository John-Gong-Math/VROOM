#!/usr/bin/env bash
# Benchmark G1 operations: VROOM (C++) vs midnight-zk (Rust/BLST) on the same GCP VM.
# Prerequisites: gcloud CLI installed and authenticated.
#
# Usage: ./bench_g1_comparison_gcp.sh [--keep]
#   --keep   Don't delete the VM after benchmarking (for debugging)

set -euo pipefail

# ── Configuration ──────────────────────────────────────────────────────────
VM_NAME="vroom-bench"
ZONE="us-central1-a"
MACHINE_TYPE="c3-standard-4"
IMAGE_FAMILY="debian-12"
IMAGE_PROJECT="debian-cloud"
KEEP_VM=false

VROOM_DIR="$(cd "$(dirname "$0")" && pwd)"
MIDNIGHT_DIR="$VROOM_DIR/../../rust/midnight-zk"

if [[ "${1:-}" == "--keep" ]]; then
    KEEP_VM=true
fi

# ── Preflight checks ──────────────────────────────────────────────────────
if ! command -v gcloud &>/dev/null; then
    echo "ERROR: gcloud CLI not found. Install it first:"
    echo "  brew install --cask google-cloud-sdk"
    exit 1
fi

if [[ ! -d "$MIDNIGHT_DIR/curves" ]]; then
    echo "ERROR: midnight-zk repo not found at $MIDNIGHT_DIR"
    echo "Expected: $MIDNIGHT_DIR/curves/Cargo.toml"
    exit 1
fi

PROJECT=$(gcloud config get-value project 2>/dev/null)
if [[ -z "$PROJECT" ]]; then
    echo "ERROR: No GCP project set. Run: gcloud init"
    exit 1
fi
echo "Using GCP project: $PROJECT"

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
    --boot-disk-size=30GB \
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

# ── Install dependencies ─────────────────────────────────────────────────
echo "Installing build dependencies..."
gcloud compute ssh "$VM_NAME" --zone="$ZONE" --command="
set -e
sudo apt-get update -qq
sudo apt-get install -y -qq clang llvm lld g++ libgmp-dev libbenchmark-dev make curl build-essential
# Install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain 1.90.0
"

# ── Upload source code ────────────────────────────────────────────────────
echo "Uploading VROOM source code..."
tar czf /tmp/vroom-src.tar.gz \
    -C "$VROOM_DIR" \
    --exclude='.git' \
    --exclude='*.o' \
    --exclude='*.a' \
    --exclude='notes' \
    .
gcloud compute scp /tmp/vroom-src.tar.gz "$VM_NAME":~/vroom-src.tar.gz --zone="$ZONE" --quiet
rm /tmp/vroom-src.tar.gz

echo "Uploading midnight-zk source code..."
tar czf /tmp/midnight-src.tar.gz \
    -C "$MIDNIGHT_DIR" \
    --exclude='.git' \
    --exclude='target' \
    .
gcloud compute scp /tmp/midnight-src.tar.gz "$VM_NAME":~/midnight-src.tar.gz --zone="$ZONE" --quiet
rm /tmp/midnight-src.tar.gz

# ── Run VROOM benchmark ──────────────────────────────────────────────────
echo ""
echo "=========================================="
echo "  Part 1: VROOM (C++ / RNS AVX512-IFMA)"
echo "=========================================="
gcloud compute ssh "$VM_NAME" --zone="$ZONE" --command="
set -e

if grep -q avx512ifma /proc/cpuinfo; then
    echo 'AVX512-IFMA: supported'
else
    echo 'WARNING: AVX512-IFMA not detected!'
    grep -m1 'model name' /proc/cpuinfo
fi

mkdir -p ~/vroom && cd ~/vroom
tar xzf ~/vroom-src.tar.gz

echo ''
echo '=== Building BLST ==='
cd ~/vroom/blst && make CC=clang

echo ''
echo '=== Building VROOM benchmark ==='
cd ~/vroom/src && make CXX=clang++ bench_pairing_50bit

echo ''
echo '=========================================='
echo '  VROOM G1 Benchmark Results'
echo '=========================================='
~/vroom/src/bench_pairing_50bit --benchmark_filter='BM_G1_'
"

# ── Run midnight-zk benchmark ────────────────────────────────────────────
echo ""
echo "=========================================="
echo "  Part 2: midnight-zk (Rust / BLST)"
echo "=========================================="
gcloud compute ssh "$VM_NAME" --zone="$ZONE" --command='
set -e
source "$HOME/.cargo/env"

mkdir -p ~/midnight && cd ~/midnight
tar xzf ~/midnight-src.tar.gz

echo ""
echo "=== Building and running midnight-zk ec benchmark ==="
echo "(This compiles with opt-level=3, LTO, codegen-units=1)"
echo ""

cd ~/midnight
cargo bench --bench ec -- "Bn256-G1 scalar multiplication" --sample-size 1000 --warm-up-time 3 --measurement-time 10
echo ""
echo "=== Running remaining G1 benchmarks ==="
cargo bench --bench ec -- "Bn256-G1 addition" --sample-size 1000 --warm-up-time 3 --measurement-time 5
cargo bench --bench ec -- "Bn256-G1 mixed addition" --sample-size 1000 --warm-up-time 3 --measurement-time 5
cargo bench --bench ec -- "Bn256-G1 doubling" --sample-size 1000 --warm-up-time 3 --measurement-time 5
'

echo ""
echo "=========================================="
echo "  Comparison Benchmark Complete"
echo "=========================================="
echo "Both benchmarks ran on the same $MACHINE_TYPE VM in $ZONE."
