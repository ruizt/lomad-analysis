#!/bin/bash
# submit_sweep.sh — submit one Kubernetes job per validation experiment
#
# Usage: bash simulations/validation/tide/submit_sweep.sh
#
# Three jobs total (see design.md for details):
#   rho-s100   (Panel A: rho accuracy)
#   var-s100   (Panel A: V accuracy)
#   e2e-s100   (Panel B: end-to-end pipeline)

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

# Panel A: Proposition 1 moment accuracy
# rho runs at 100 reps, not 2000: at 500 the Monte Carlo error is smaller than
# the line width and the empirical curve disappears under the theoretical one.
submit_job "rho-s100"  100   1  "1Gi"
submit_job "var-s100"  2000  1  "2Gi"

# Panel B: end-to-end pipeline
submit_job "e2e-s100"  1000  2  "2Gi"

echo ""
echo "All jobs submitted. Monitor with:"
echo "  kubectl get jobs -n ${NAMESPACE} -l app=lomad-valid"
echo "  kubectl logs -n ${NAMESPACE} job/<job-name>"
