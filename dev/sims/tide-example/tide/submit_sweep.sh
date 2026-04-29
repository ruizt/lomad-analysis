#!/bin/bash
# submit_sweep.sh — submit one Kubernetes job per sample size n
#
# Usage: bash dev/sims/tide-example/tide/submit_sweep.sh
#
# Each job runs SIM_S replicates for one value of n.
# Results land in the mvn-example-results PVC as one .rds file per job
# (e.g. n10.rds, n30.rds, n100.rds, n500.rds).
# Collect after all jobs complete with tide/collect.R.

NAMESPACE="cal-poly-lomad"
N_VALUES=(10 30 100 500)
SIM_S=200
SIM_SEED=7291
IMAGE="ghcr.io/ruizt/lomad-sims:latest"
CONFIGMAP="mvn-example-script"

# Create/update the ConfigMap from the local sim.R
echo "Creating ConfigMap '${CONFIGMAP}' ..."
kubectl create configmap ${CONFIGMAP} \
  -n ${NAMESPACE} \
  --from-file=sim.R=dev/sims/tide-example/tide/sim.R \
  --dry-run=client -o yaml | kubectl apply -f -
echo ""

for n in "${N_VALUES[@]}"; do
  job_name="mvn-example-n${n}"

  echo "Submitting ${job_name} (n=${n}) ..."

  kubectl apply -n ${NAMESPACE} -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: ${job_name}
  namespace: ${NAMESPACE}
  labels:
    app: mvn-example
    n: "${n}"
spec:
  backoffLimit: 1
  template:
    spec:
      restartPolicy: Never
      containers:
        - name: mvn-example
          image: ${IMAGE}
          imagePullPolicy: Always
          resources:
            requests:
              cpu: "500m"
              memory: "512Mi"
            limits:
              cpu: "500m"
              memory: "512Mi"
          env:
            - name: SIM_N
              value: "${n}"
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
            claimName: mvn-example-results
EOF

done

echo ""
echo "All jobs submitted. Monitor with:"
echo "  kubectl get jobs -n ${NAMESPACE} -l app=mvn-example"
echo "  kubectl logs -n ${NAMESPACE} job/<job-name>"
