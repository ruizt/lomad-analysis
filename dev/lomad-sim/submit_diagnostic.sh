#!/bin/bash
# submit_diagnostic.sh — two parallel sweeps to diagnose CLT rejection behaviour
#
# Sweep A (auto-h): extreme d values with auto window selection.
#   Tests whether the CLT ever rejects as d grows very large.
#
# Sweep B (h=30):   moderate d values with a fixed h=30 (matching the old
#   bootstrap setting) to test whether h=5 is hurting power.
#
# All jobs write to the same shared PVC and append to the same CSV log,
# so results from both sweeps are directly comparable.
#
# Usage:
#   bash dev/lomad-sim/submit_diagnostic.sh

set -euo pipefail

# Canonical image — built and pushed automatically from the main branch.
IMAGE="ghcr.io/ruizt/lomad-sim:latest"
NAMESPACE="cal-poly-lomad"

echo "Using image: ${IMAGE}"
echo ""

# ── helper: submit one job ────────────────────────────────────────────────────
submit_job() {
  local job_name="$1"
  local sim_d="$2"
  local sim_h="$3"
  local sim_s="$4"

  echo "  Submitting ${job_name}  (d=${sim_d}, h=${sim_h}, s=${sim_s})..."

  kubectl apply -f - --namespace="${NAMESPACE}" <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: ${job_name}
  labels:
    app: lomad-diag
spec:
  backoffLimit: 1
  ttlSecondsAfterFinished: 7200
  template:
    spec:
      restartPolicy: Never
      containers:
        - name: lomad-sim
          image: ${IMAGE}
          imagePullPolicy: Always
          resources:
            requests:
              cpu: "1"
              memory: "1Gi"
            limits:
              cpu: "1200m"
              memory: "1228Mi"
          env:
            - name: SIM_N
              value: "500"
            - name: SIM_D
              value: "${sim_d}"
            - name: SIM_SEED
              value: "32026"
            - name: SIM_H
              value: "${sim_h}"
            - name: SIM_S
              value: "${sim_s}"
            - name: SIM_MAX_PQ
              value: "3"
            - name: SIM_OUT_DIR
              value: "/jobs/output"
          volumeMounts:
            - name: output
              mountPath: /jobs/output
      volumes:
        - name: output
          persistentVolumeClaim:
            claimName: lomad-sim-results
EOF
}

# ── Sweep A: extreme d, h=auto ────────────────────────────────────────────────
# Tests whether the CLT ever rejects as separation grows very large.
# h=0 and s=0 trigger auto-selection inside sim.R.
echo "Sweep A — extreme d, h=auto:"
for d in 10 25 50 100; do
  job_name="lomad-diag-auto-d${d}"
  submit_job "${job_name}" "${d}" "0" "0"
done

echo ""

# ── Sweep B: moderate d, h=30 ─────────────────────────────────────────────────
# Tests whether h=5 (auto) was the reason for zero rejections.
# h=30 matches the window the old bootstrap method used.
# s is set to 0 (auto) so only h changes relative to Sweep A.
echo "Sweep B — moderate d, h=30:"
for d in 5 10 25 50; do
  job_name="lomad-diag-h30-d${d}"
  submit_job "${job_name}" "${d}" "30" "0"
done

echo ""
echo "All diagnostic jobs submitted ($(( 4 + 4 )) total)."
echo "Monitor with:"
echo "  kubectl get jobs -l app=lomad-diag -w --namespace=${NAMESPACE}"
echo "  kubectl get pods -w --namespace=${NAMESPACE}"
echo ""
echo "When complete, run:"
echo "  bash lomad-sim/fetch_results.sh"
