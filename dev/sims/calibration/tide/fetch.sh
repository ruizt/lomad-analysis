#!/bin/bash
# fetch.sh — copy .rds results from the PVC to a local directory
#
# Usage (from the repo root):
#   bash dev/sims/calibration/tide/fetch.sh
#
# Optional: override the local destination as the first argument:
#   bash dev/sims/calibration/tide/fetch.sh /tmp/my-results
#
# How it works:
#   1. Spins up a lightweight accessor pod that mounts the lomad-calib-results PVC.
#   2. Waits for the pod to be Ready (up to 2 minutes).
#   3. kubectl cp copies all .rds files to LOCAL_DIR.
#   4. Deletes the accessor pod.
#
# After this completes, assemble the results in R:
#   Rscript dev/sims/calibration/tide/collect.R

set -euo pipefail

LOCAL_DIR="${1:-dev/sims/calibration/results/raw}"
ACCESSOR_YAML="dev/sims/calibration/tide/accessor.yaml"
POD_NAME="lomad-fetch"
NAMESPACE="cal-poly-lomad"

echo "=== lomad calibration: fetch results from PVC ==="
echo "Destination: ${LOCAL_DIR}"
echo ""

# Create local directory
mkdir -p "${LOCAL_DIR}"

# Start accessor pod
echo "[1/4] Starting accessor pod ..."
kubectl apply -n "${NAMESPACE}" -f "${ACCESSOR_YAML}"

# Wait until the pod is Running
echo "[2/4] Waiting for pod to be ready ..."
kubectl wait -n "${NAMESPACE}" --for=condition=Ready "pod/${POD_NAME}" --timeout=120s

# Copy all .rds files
echo "[3/4] Copying .rds files ..."
kubectl cp -n "${NAMESPACE}" "${POD_NAME}:/jobs/output/." "${LOCAL_DIR}/"
echo "      Done."

# Tear down accessor pod
echo "[4/4] Cleaning up accessor pod ..."
kubectl delete -n "${NAMESPACE}" pod "${POD_NAME}" --ignore-not-found

echo ""
echo "Results are in: ${LOCAL_DIR}"
echo ""
echo "Next step — assemble into summary.rds:"
echo "  Rscript dev/sims/calibration/tide/collect.R"
