#!/bin/bash
# test_one_job.sh — submit a single job to verify sim.R works on the cluster
#
# Usage: bash simulations/power/tide/test_one_job.sh
#
# Uses the same tide/job.yaml template as submit_sweep.sh, so a job submitted
# here is identical to a sweep job apart from the parameters set below. Edit
# those to probe a different cell; SIM_S is small so it finishes quickly.

set -euo pipefail

command -v envsubst >/dev/null 2>&1 || {
  echo "error: envsubst not found (part of gettext)." >&2
  echo "       install with: brew install gettext" >&2
  exit 1
}

TIDE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_TEMPLATE="${TIDE_DIR}/job.yaml"

export NAMESPACE="cal-poly-ruiz"
export IMAGE="ghcr.io/ruizt/lomad-simulations:latest"
export CONFIGMAP="lomad-power-script"

# ---- The cell to test --------------------------------------------------------

export SIM_STRUCTURE="rate"
export SIM_D="1.0"
export SIM_N="400"
export SIM_SNR="0.5"
export SIM_PHI="0.8"
export SIM_ORACLE="FALSE"
export SIM_S=5          # deliberately small: this is a smoke test
export SIM_SEED=2847

export JOB_NAME="lomad-power-test"

# ---- Submit ------------------------------------------------------------------

echo "Updating ConfigMap '${CONFIGMAP}' ..."
kubectl create configmap "${CONFIGMAP}" \
  -n "${NAMESPACE}" \
  --from-file=sim.R="${TIDE_DIR}/sim.R" \
  --dry-run=client -o yaml | kubectl apply -f -

echo "Submitting ${JOB_NAME} (S=${SIM_S}) ..."
envsubst '${JOB_NAME} ${NAMESPACE} ${IMAGE} ${CONFIGMAP} ${SIM_D} ${SIM_STRUCTURE} ${SIM_N} ${SIM_SNR} ${SIM_PHI} ${SIM_S} ${SIM_SEED} ${SIM_ORACLE}' \
  < "${JOB_TEMPLATE}" | kubectl apply -n "${NAMESPACE}" -f -

echo ""
echo "Watching job status (Ctrl-C to stop) ..."
kubectl get job "${JOB_NAME}" -n "${NAMESPACE}" -w

echo ""
echo "Logs:    kubectl logs -n ${NAMESPACE} job/${JOB_NAME}"
echo "Clean up: kubectl delete job -n ${NAMESPACE} ${JOB_NAME}"
