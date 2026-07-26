#!/bin/bash
# fetch.sh — copy .rds results from the PVC to a local directory
#
# Usage (from the repo root):
#   bash simulations/validation/tide/fetch.sh
#
# Optional: override the local destination as the first argument:
#   bash simulations/validation/tide/fetch.sh /tmp/my-results
#
# After this completes, assemble the results in R:
#   Rscript simulations/validation/collect.R

set -euo pipefail

LOCAL_DIR="${1:-simulations/validation/results/raw}"
ACCESSOR_YAML="simulations/validation/tide/accessor.yaml"
POD_NAME="lomad-valid-fetch"
NAMESPACE="cal-poly-ruiz"

echo "=== lomad validation: fetch results from PVC ==="
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
echo "Next step — assemble results and generate figures:"
echo "  Rscript simulations/validation/collect.R"
