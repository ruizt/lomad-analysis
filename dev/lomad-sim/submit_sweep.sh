#!/bin/bash
# submit_sweep.sh — submit one Kubernetes job per value of d
#
# Usage:
#   bash lomad-sim/submit_sweep.sh
#
# To use a different GitHub username (e.g. a collaborator's fork):
#   GHCR_USER=their-username bash lomad-sim/submit_sweep.sh

# ── CONFIGURE ────────────────────────────────────────────────────────────────
# Your GitHub username — controls which container registry image is pulled.
# Override at runtime: GHCR_USER=yourname bash submit_sweep.sh
GHCR_USER="${GHCR_USER:-otishunt}"

# Separation values to sweep over
D_VALUES=(0.5 1 2 3 5)
# ─────────────────────────────────────────────────────────────────────────────

IMAGE="ghcr.io/${GHCR_USER}/lomad-sim:latest"
echo "Using image: ${IMAGE}"

for d in "${D_VALUES[@]}"; do
  # Kubernetes job names must be DNS-safe (no dots)
  job_name="lomad-sim-d$(echo $d | tr '.' '-')"

  echo "Submitting $job_name (d=$d)..."

  kubectl apply -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: $job_name
  labels:
    app: lomad-sim
    d: "$d"
spec:
  backoffLimit: 1
  # ttlSecondsAfterFinished keeps the pod (and its logs) alive for 2 hours
  # after completion so you can inspect them before they are cleaned up.
  ttlSecondsAfterFinished: 7200
  template:
    spec:
      restartPolicy: Never
      imagePullSecrets:
        - name: ghcr-secret
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
              value: "$d"
            - name: SIM_SEED
              value: "32026"
            - name: SIM_H
              value: "0"
            - name: SIM_S
              value: "0"
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

done

echo ""
echo "All jobs submitted. Monitor with:"
echo "  kubectl get jobs -w"
echo "  kubectl get pods -w"
