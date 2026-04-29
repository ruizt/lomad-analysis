#!/bin/bash
# fetch.sh — copy .rds results from the PVC to a local directory
#
# Usage (from the repo root):
#   bash dev/sims/power/tide/fetch.sh
#
# Optional: override the local destination as the first argument:
#   bash dev/sims/power/tide/fetch.sh /tmp/my-results

set -euo pipefail

LOCAL_DIR="${1:-dev/sims/power/results/raw}"
ACCESSOR_YAML="dev/sims/power/tide/accessor.yaml"
POD_NAME="power-fetch"
NAMESPACE="cal-poly-lomad"

echo "=== lomad power: fetch results from PVC ==="
echo "Destination: ${LOCAL_DIR}"
echo ""

mkdir -p "${LOCAL_DIR}"

echo "[1/4] Starting accessor pod ..."
kubectl apply -n "${NAMESPACE}" -f "${ACCESSOR_YAML}"

echo "[2/4] Waiting for pod to be ready ..."
kubectl wait -n "${NAMESPACE}" --for=condition=Ready "pod/${POD_NAME}" --timeout=120s

echo "[3/4] Copying .rds files ..."
kubectl cp -n "${NAMESPACE}" "${POD_NAME}:/jobs/output/." "${LOCAL_DIR}/"
echo "      Done."

echo "[4/4] Cleaning up accessor pod ..."
kubectl delete -n "${NAMESPACE}" pod "${POD_NAME}" --ignore-not-found

echo ""
echo "Results are in: ${LOCAL_DIR}"
echo ""
echo "Next step — assemble results:"
echo "  Rscript dev/sims/power/tide/collect.R"
