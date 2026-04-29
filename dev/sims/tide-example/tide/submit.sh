#!/bin/bash
# submit.sh — one-shot script to run the full MVN example on Tide
#
# Usage (from the repo root):
#   bash dev/sims/tide-example/tide/submit.sh
#
# What it does:
#   1. Creates the PVC (idempotent — safe to re-run)
#   2. Submits the sweep (creates ConfigMap + Jobs)
#   3. Waits for all jobs to complete
#   4. Fetches results from the PVC
#   5. Cleans up jobs and accessor pod
#
# The PVC is left in place so you can inspect results later.
# To delete it: kubectl delete -n cal-poly-lomad -f dev/sims/tide-example/tide/pvc.yaml

set -euo pipefail

NAMESPACE="cal-poly-lomad"
STUDY_DIR="dev/sims/tide-example"

echo "=== MVN coverage example: full pipeline ==="
echo ""

# ---- Step 1: Create PVC -----------------------------------------------------

echo "[1/4] Creating PVC ..."
kubectl apply -n ${NAMESPACE} -f ${STUDY_DIR}/tide/pvc.yaml
echo ""

# ---- Step 2: Submit jobs -----------------------------------------------------

echo "[2/4] Submitting jobs ..."
bash ${STUDY_DIR}/tide/submit_sweep.sh
echo ""

# ---- Step 3: Wait for completion ---------------------------------------------

echo "[3/4] Waiting for jobs to complete ..."
while true; do
  total=$(kubectl get jobs -n ${NAMESPACE} -l app=mvn-example --no-headers 2>/dev/null | wc -l | tr -d ' ')
  complete=$(kubectl get jobs -n ${NAMESPACE} -l app=mvn-example --no-headers 2>/dev/null | grep -c "1/1" || true)

  echo "      ${complete}/${total} complete"

  if [ "${complete}" -eq "${total}" ] && [ "${total}" -gt 0 ]; then
    break
  fi

  sleep 15
done
echo ""

# ---- Step 4: Fetch results ---------------------------------------------------

echo "[4/4] Fetching results ..."
bash ${STUDY_DIR}/tide/fetch.sh
echo ""

# ---- Done --------------------------------------------------------------------

echo "=== Pipeline complete ==="
echo ""
echo "Assemble results locally with:"
echo "  Rscript ${STUDY_DIR}/tide/collect.R"
echo ""
echo "Clean up cluster resources with:"
echo "  kubectl delete jobs -n ${NAMESPACE} -l app=mvn-example"
echo "  kubectl delete -n ${NAMESPACE} -f ${STUDY_DIR}/tide/pvc.yaml"
