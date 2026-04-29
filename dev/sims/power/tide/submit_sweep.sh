#!/bin/bash
# submit_sweep.sh — submit one Kubernetes job per (structure, d) combination
#
# Usage: bash dev/sims/power/tide/submit_sweep.sh
#
# Each job runs SIM_S replicates for one (structure, d) pair.
# Results land in the lomad-power-results PVC as one .rds file per job
# (e.g. dist_d0-0.rds, smooth_d1-5.rds).
# Collect after all jobs complete with tide/collect.R.

NAMESPACE="cal-poly-lomad"
D_VALUES=(0 0.5 1.0 1.5 2.0)
STRUCTURES=(dist smooth cross rate)
SIM_S=200
SIM_SEED=2847
IMAGE="ghcr.io/ruizt/lomad-sims:latest"
CONFIGMAP="lomad-power-script"

# Create/update the ConfigMap from the local sim.R
echo "Creating ConfigMap '${CONFIGMAP}' ..."
kubectl create configmap ${CONFIGMAP} \
  -n ${NAMESPACE} \
  --from-file=sim.R=dev/sims/power/tide/sim.R \
  --dry-run=client -o yaml | kubectl apply -f -
echo ""

for structure in "${STRUCTURES[@]}"; do
  for d in "${D_VALUES[@]}"; do
    d_label=$(echo $d | tr '.' '-')
    job_name="lomad-power-${structure}-d${d_label}"

    echo "Submitting ${job_name} (structure=${structure}, d=${d}) ..."

    kubectl apply -n ${NAMESPACE} -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: ${job_name}
  namespace: ${NAMESPACE}
  labels:
    app: lomad-power
    structure: "${structure}"
    d: "${d}"
spec:
  backoffLimit: 1
  template:
    spec:
      restartPolicy: Never
      containers:
        - name: lomad-power
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
            - name: SIM_STRUCTURE
              value: "${structure}"
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
            claimName: lomad-power-results
EOF

  done
done

echo ""
echo "All jobs submitted. Monitor with:"
echo "  kubectl get jobs -n ${NAMESPACE} -l app=lomad-power"
echo "  kubectl logs -n ${NAMESPACE} job/<job-name>"
