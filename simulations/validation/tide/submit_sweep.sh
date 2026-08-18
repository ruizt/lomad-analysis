#!/bin/bash
# submit_sweep.sh — submit one Kubernetes job per validation experiment
#
# Usage: bash simulations/validation/tide/submit_sweep.sh
#
# One job (see design.md for details):
#   e2e-s100   both panels: Proposition 1 moments and pipeline calibration

NAMESPACE="cal-poly-ruiz"
SIM_SEED=7291
IMAGE="ghcr.io/ruizt/lomad-simulations:latest"
CONFIGMAP="lomad-valid-script"

# Create/update the ConfigMap from the local sim.R
echo "Creating ConfigMap '${CONFIGMAP}' ..."
kubectl create configmap ${CONFIGMAP} \
  -n ${NAMESPACE} \
  --from-file=sim.R=simulations/validation/tide/sim.R \
  --dry-run=client -o yaml | kubectl apply -f -
echo ""

# Helper: submit one job
submit_job() {
  local exp="$1"
  local sim_s="$2"
  local cpu="$3"
  local mem="$4"

  local job_name="lomad-valid-${exp}"

  echo "Submitting ${job_name} (S=${sim_s}, cpu=${cpu}, mem=${mem}) ..."

  kubectl apply -n ${NAMESPACE} -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: ${job_name}
  namespace: ${NAMESPACE}
  labels:
    app: lomad-valid
    experiment: "${exp}"
spec:
  backoffLimit: 1
  template:
    spec:
      restartPolicy: Never
      containers:
        - name: lomad-valid
          image: ${IMAGE}
          imagePullPolicy: Always
          resources:
            requests:
              cpu: "${cpu}"
              memory: "${mem}"
            limits:
              cpu: "${cpu}"
              memory: "${mem}"
          env:
            - name: SIM_EXPERIMENT
              value: "${exp}"
            - name: SIM_S
              value: "${sim_s}"
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
            claimName: lomad-valid-results
EOF
}

submit_job "e2e-s100"  1000  2  "2Gi"

echo ""
echo "All jobs submitted. Monitor with:"
echo "  kubectl get jobs -n ${NAMESPACE} -l app=lomad-valid"
echo "  kubectl logs -n ${NAMESPACE} job/<job-name>"
