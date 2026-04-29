#!/bin/bash
# submit_sweep.sh — submit one Kubernetes job per value of d
#
# Usage: bash dev/sims/calibration/tide/submit_sweep.sh
#
# Each job runs SIM_S replicates for one value of d.
# Results land in the lomad-calib-results PVC as one .rds file per job
# (e.g. d0-0.rds, d0-1.rds, ..., d0-5.rds).
# Collect after all jobs complete with tide/collect.R.

D_VALUES=(0 0.1 0.2 0.3 0.4 0.5)
SIM_S=200
SIM_SEED=4853
IMAGE="ghcr.io/<org>/lomad-calib:latest"

for d in "${D_VALUES[@]}"; do
  job_name="lomad-calib-d$(echo $d | tr '.' '-')"

  echo "Submitting ${job_name} (d=${d}) ..."

  kubectl apply -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: ${job_name}
  labels:
    app: lomad-calib
    d: "${d}"
spec:
  backoffLimit: 1
  template:
    spec:
      restartPolicy: Never
      imagePullSecrets:
        - name: ghcr-secret
      containers:
        - name: lomad-calib
          image: ${IMAGE}
          imagePullPolicy: Always
          resources:
            requests:
              cpu: "1"
              memory: "2Gi"
            limits:
              cpu: "1200m"
              memory: "2Gi"
          env:
            - name: SIM_D
              value: "${d}"
            - name: SIM_S
              value: "${SIM_S}"
            - name: SIM_SEED
              value: "${SIM_SEED}"
            - name: SIM_OUT_DIR
              value: "/jobs/output"
          volumeMounts:
            - name: output
              mountPath: /jobs/output
      volumes:
        - name: output
          persistentVolumeClaim:
            claimName: lomad-calib-results
EOF

done

echo ""
echo "All jobs submitted. Monitor with:"
echo "  kubectl get jobs -l app=lomad-calib"
echo "  kubectl logs job/<job-name>"
