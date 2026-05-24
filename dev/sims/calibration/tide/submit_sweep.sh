#!/bin/bash
# submit_sweep.sh — submit one Kubernetes job per value of d
#
# Usage: bash dev/sims/calibration/tide/submit_sweep.sh
#
# Each job runs SIM_S replicates for one value of d.
# Results land in the lomad-calib-results PVC as one .rds file per job
# (e.g. d0-0.rds, d0-2.rds, d0-5.rds, d1-0.rds).
# Collect after all jobs complete with tide/collect.R.

NAMESPACE="cal-poly-ruiz"
D_VALUES=(0.25 0.5 0.75 1 1.25 1.5 1.75 2)
SIM_S=500
SIM_SEED=4853
IMAGE="ghcr.io/ruizt/lomad-sims:latest"
CONFIGMAP="lomad-calib-script"

# Create/update the ConfigMap from the local sim.R
echo "Creating ConfigMap '${CONFIGMAP}' ..."
kubectl create configmap ${CONFIGMAP} \
  -n ${NAMESPACE} \
  --from-file=sim.R=dev/sims/calibration/tide/sim.R \
  --dry-run=client -o yaml | kubectl apply -f -
echo ""

for d in "${D_VALUES[@]}"; do
  job_name="lomad-calib-d$(echo $d | tr '.' '-')"

  echo "Submitting ${job_name} (d=${d}) ..."

  kubectl apply -n ${NAMESPACE} -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: ${job_name}
  namespace: ${NAMESPACE}
  labels:
    app: lomad-calib
    d: "${d}"
spec:
  backoffLimit: 1
  template:
    spec:
      restartPolicy: Never
      containers:
        - name: lomad-calib
          image: ${IMAGE}
          imagePullPolicy: Always
          resources:
            requests:
              cpu: "1"
              memory: "1Gi"
            limits:
              cpu: "1"
              memory: "1Gi"
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
            - name: script
              mountPath: /scripts
            - name: output
              mountPath: /jobs/output
      volumes:
        - name: script
          configMap:
            name: ${CONFIGMAP}
        - name: output
          persistentVolumeClaim:
            claimName: lomad-calib-results
EOF

done

echo ""
echo "All jobs submitted. Monitor with:"
echo "  kubectl get jobs -n ${NAMESPACE} -l app=lomad-calib"
echo "  kubectl logs -n ${NAMESPACE} job/<job-name>"
