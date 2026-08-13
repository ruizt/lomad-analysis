#!/bin/bash
# submit_sweep.sh — submit one Kubernetes Job per (structure, d, s, snr, phi)
#
# Usage: bash simulations/power/tide/submit_sweep.sh
#        STRUCTS="fr" bash simulations/power/tide/submit_sweep.sh   # subset
#
# Each Job runs SIM_REPS replicates for one parameter combination. Results land
# in the lomad-power-results PVC as one .rds file per Job (plus a -windows.rds).
# Fetch with tide/fetch.sh, then assemble with collect-results.R.
#
# The Job spec lives in tide/job.yaml and is filled in here with envsubst, so
# there is exactly one copy of it to keep in step with sim.R.

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
export SIM_REPS=500
export SIM_SEED=2847

# Structures to sweep. Override to rerun a subset in place -- results are one
# file per cell in the PVC, so untouched structures are left alone:
#   STRUCTS="fr" bash simulations/power/tide/submit_sweep.sh
read -r -a STRUCTURES <<< "${STRUCTS:-rs rm fr}"

# The window s is the design factor; sim.R derives n = 25 s from it. Holding
# n/s fixed keeps the accumulated affine drift identical across window sizes
# instead of confounding the two.
WINDOW_SIZES=(50 100 150)
SNR_VALUES=(0.5 1.5)
PHI_VALUES=(0.3 0.5 0.7)

# d is a generator knob, not an effect size: the same d gives ~2x different
# local separation across structures. It is scaled per structure onto a common
# realized delta_t (median ~0, 0.20, 0.45). See _notes/delta-calibration.md.
d_values_for() {
  case "$1" in
    rs) echo "0 0.60 1.55" ;;
    rm) echo "0 0.36 0.93" ;;
    fr) echo "0 0.5 1.7" ;;
    *)  echo "error: no d grid for structure '$1'" >&2; exit 1 ;;
  esac
}

# ---- ConfigMap ---------------------------------------------------------------

echo "Creating ConfigMap '${CONFIGMAP}' ..."
kubectl create configmap "${CONFIGMAP}" \
  -n "${NAMESPACE}" \
  --from-file=sim.R="${TIDE_DIR}/sim.R" \
  --dry-run=client -o yaml | kubectl apply -f -
echo ""

# ---- Submit one Job from the template ----------------------------------------
# Expects SIM_STRUCTURE, SIM_D, SIM_S, SIM_SNR, SIM_PHI exported.

submit_job() {
  local label_d label_snr label_phi
  label_d=${SIM_D//./-}
  label_snr=${SIM_SNR//./-}
  label_phi=${SIM_PHI//./-}

  export JOB_NAME="lomad-power-${SIM_STRUCTURE}-d${label_d}-s${SIM_S}-snr${label_snr}-phi${label_phi}"

  echo "Submitting ${JOB_NAME} ..."
  envsubst '${JOB_NAME} ${NAMESPACE} ${IMAGE} ${CONFIGMAP} ${SIM_D} ${SIM_STRUCTURE} ${SIM_S} ${SIM_SNR} ${SIM_PHI} ${SIM_REPS} ${SIM_SEED}' \
    < "${JOB_TEMPLATE}" | kubectl apply -n "${NAMESPACE}" -f -
}

# ---- Main sweep --------------------------------------------------------------

for SIM_STRUCTURE in "${STRUCTURES[@]}"; do
  read -r -a D_VALUES <<< "$(d_values_for "${SIM_STRUCTURE}")"
  for SIM_D in "${D_VALUES[@]}"; do
    for SIM_S in "${WINDOW_SIZES[@]}"; do
      for SIM_SNR in "${SNR_VALUES[@]}"; do
        for SIM_PHI in "${PHI_VALUES[@]}"; do
          export SIM_STRUCTURE SIM_D SIM_S SIM_SNR SIM_PHI
          submit_job
        done
      done
    done
  done
done

echo ""
echo "All jobs submitted. Monitor with:"
echo "  kubectl get jobs -n ${NAMESPACE} -l app=lomad-power"
echo "  kubectl logs -n ${NAMESPACE} job/<job-name>"
