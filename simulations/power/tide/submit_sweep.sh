#!/bin/bash
# submit_sweep.sh — submit one Kubernetes Job per (structure, d, n, snr, phi)
#
# Usage: bash simulations/power/tide/submit_sweep.sh
#
# Each Job runs SIM_S replicates for one parameter combination. Results land in
# the lomad-power-results PVC as one .rds file per Job (plus a -series.rds).
# Fetch with tide/fetch.sh, then assemble with collect.R.
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
export IMAGE="ghcr.io/ruizt/lomad-sims:latest"
export CONFIGMAP="lomad-power-script"
export SIM_S=500
export SIM_SEED=2847

D_VALUES=(0 0.5 1.0 1.5 2.0)
STRUCTURES=(smooth cross rate)
SAMPLE_SIZES=(200 400 600)
SNR_VALUES=(0.5 1.5)
PHI_VALUES=(0.3 0.5 0.8)

# ---- ConfigMap ---------------------------------------------------------------

echo "Creating ConfigMap '${CONFIGMAP}' ..."
kubectl create configmap "${CONFIGMAP}" \
  -n "${NAMESPACE}" \
  --from-file=sim.R="${TIDE_DIR}/sim.R" \
  --dry-run=client -o yaml | kubectl apply -f -
echo ""

# ---- Submit one Job from the template ----------------------------------------
# Expects SIM_STRUCTURE, SIM_D, SIM_N, SIM_SNR, SIM_PHI, SIM_ORACLE exported.

submit_job() {
  local label_d label_snr label_phi suffix
  label_d=${SIM_D//./-}
  label_snr=${SIM_SNR//./-}
  label_phi=${SIM_PHI//./-}
  suffix=""
  [ "${SIM_ORACLE}" = "TRUE" ] && suffix="-oracle"

  export JOB_NAME="lomad-power-${SIM_STRUCTURE}-d${label_d}-n${SIM_N}-snr${label_snr}-phi${label_phi}${suffix}"

  echo "Submitting ${JOB_NAME} ..."
  envsubst '${JOB_NAME} ${NAMESPACE} ${IMAGE} ${CONFIGMAP} ${SIM_D} ${SIM_STRUCTURE} ${SIM_N} ${SIM_SNR} ${SIM_PHI} ${SIM_S} ${SIM_SEED} ${SIM_ORACLE}' \
    < "${JOB_TEMPLATE}" | kubectl apply -n "${NAMESPACE}" -f -
}

# ---- Main sweep (estimated noise) --------------------------------------------

export SIM_ORACLE="FALSE"

for SIM_STRUCTURE in "${STRUCTURES[@]}"; do
  for SIM_D in "${D_VALUES[@]}"; do
    for SIM_N in "${SAMPLE_SIZES[@]}"; do
      for SIM_SNR in "${SNR_VALUES[@]}"; do
        for SIM_PHI in "${PHI_VALUES[@]}"; do
          export SIM_STRUCTURE SIM_D SIM_N SIM_SNR SIM_PHI
          submit_job
        done
      done
    done
  done
done

# ---- Oracle sweep (phi = 0.8 only) -------------------------------------------
# Isolates estimation error from test behaviour: sim.R bypasses noise
# estimation and uses the true AR(1) parameters. See design.md.

echo ""
echo "Submitting oracle jobs (phi = 0.8 only) ..."

export SIM_ORACLE="TRUE"
export SIM_PHI="0.8"

for SIM_STRUCTURE in "${STRUCTURES[@]}"; do
  for SIM_D in "${D_VALUES[@]}"; do
    for SIM_N in "${SAMPLE_SIZES[@]}"; do
      for SIM_SNR in "${SNR_VALUES[@]}"; do
        export SIM_STRUCTURE SIM_D SIM_N SIM_SNR
        submit_job
      done
    done
  done
done

echo ""
echo "All jobs submitted. Monitor with:"
echo "  kubectl get jobs -n ${NAMESPACE} -l app=lomad-power"
echo "  kubectl get jobs -n ${NAMESPACE} -l oracle=TRUE      # oracle arm only"
echo "  kubectl logs -n ${NAMESPACE} job/<job-name>"
