#!/bin/bash
# submit_sweep.sh — submit one Kubernetes job per validation experiment
#
# Usage: bash dev/sims/validation/tide/submit_sweep.sh
#
# Six jobs total (see design.md for details):
#   clt-s80, clt-s150, clt-s300   (Figure 1: CLT QQ)
#   rho-s150                      (Figure 2: rho accuracy)
#   var-s150                      (Figure 2: V accuracy)
#   e2e-s150                      (Figure 3: end-to-end pipeline)

NAMESPACE="cal-poly-ruiz"
SIM_SEED=7291
IMAGE="ghcr.io/ruizt/lomad-sims:latest"
CONFIGMAP="lomad-valid-script"

# Create/update the ConfigMap from the local sim.R
echo "Creating ConfigMap '${CONFIGMAP}' ..."
kubectl create configmap ${CONFIGMAP} \
  -n ${NAMESPACE} \
  --from-file=sim.R=dev/sims/validation/tide/sim.R \
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

# Figure 1: CLT QQ (3 jobs)
submit_job "clt-s80"   2000  1  "1Gi"
submit_job "clt-s150"  2000  1  "1Gi"
submit_job "clt-s300"  2000  1  "1Gi"

# Figure 2: Proposition 1 moment accuracy (2 jobs, both s = 150)
submit_job "rho-s150"  500   1  "1Gi"
submit_job "var-s150"  2000  1  "2Gi"

# Figure 3: end-to-end pipeline (1 job)
submit_job "e2e-s150"  1000  2  "2Gi"

echo ""
echo "All jobs submitted. Monitor with:"
echo "  kubectl get jobs -n ${NAMESPACE} -l app=lomad-valid"
echo "  kubectl logs -n ${NAMESPACE} job/<job-name>"
