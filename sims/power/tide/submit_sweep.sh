#!/bin/bash
# submit_sweep.sh — submit one Kubernetes job per (structure, d, n, snr) combination
#
# Usage: bash sims/power/tide/submit_sweep.sh
#
# Each job runs SIM_S replicates for one parameter combination.
# Results land in the lomad-power-results PVC as one .rds file per job
# (e.g. smooth_d0-5_n500_snr1-5.rds).
# Collect after all jobs complete with tide/collect.R.

NAMESPACE="cal-poly-ruiz"
D_VALUES=(0 0.5 1.0 1.5 2.0)
STRUCTURES=(smooth cross rate)
SAMPLE_SIZES=(200 400 600)
SNR_VALUES=(0.5 1.5)
PHI_VALUES=(0.3 0.5 0.8)
SIM_S=500
SIM_SEED=2847
IMAGE="ghcr.io/ruizt/lomad-sims:latest"
CONFIGMAP="lomad-power-script"

# Create/update the ConfigMap from the local sim.R
echo "Creating ConfigMap '${CONFIGMAP}' ..."
kubectl create configmap ${CONFIGMAP} \
  -n ${NAMESPACE} \
  --from-file=sim.R=sims/power/tide/sim.R \
  --dry-run=client -o yaml | kubectl apply -f -
echo ""

for structure in "${STRUCTURES[@]}"; do
  for d in "${D_VALUES[@]}"; do
    for n in "${SAMPLE_SIZES[@]}"; do
      for snr in "${SNR_VALUES[@]}"; do
        for phi in "${PHI_VALUES[@]}"; do
          d_label=$(echo $d | tr '.' '-')
          snr_label=$(echo $snr | tr '.' '-')
          phi_label=$(echo $phi | tr '.' '-')
          job_name="lomad-power-${structure}-d${d_label}-n${n}-snr${snr_label}-phi${phi_label}"

          echo "Submitting ${job_name} ..."

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
    n: "${n}"
    snr: "${snr}"
    phi: "${phi}"
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
            - name: SIM_N
              value: "${n}"
            - name: SIM_SNR
              value: "${snr}"
            - name: SIM_PHI
              value: "${phi}"
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
    done
  done
done

# ---- Oracle sweep (phi = 0.8 only) -------------------------------------------

echo ""
echo "Submitting oracle jobs (phi = 0.8 only) ..."

for structure in "${STRUCTURES[@]}"; do
  for d in "${D_VALUES[@]}"; do
    for n in "${SAMPLE_SIZES[@]}"; do
      for snr in "${SNR_VALUES[@]}"; do
        d_label=$(echo $d | tr '.' '-')
        snr_label=$(echo $snr | tr '.' '-')
        job_name="lomad-power-${structure}-d${d_label}-n${n}-snr${snr_label}-phi0-8-oracle"

        echo "Submitting ${job_name} ..."

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
    n: "${n}"
    snr: "${snr}"
    phi: "0.8"
    oracle: "true"
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
            - name: SIM_N
              value: "${n}"
            - name: SIM_SNR
              value: "${snr}"
            - name: SIM_PHI
              value: "0.8"
            - name: SIM_S
              value: "${SIM_S}"
            - name: SIM_SEED
              value: "${SIM_SEED}"
            - name: SIM_ORACLE
              value: "TRUE"
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
  done
done

echo ""
echo "All jobs submitted. Monitor with:"
echo "  kubectl get jobs -n ${NAMESPACE} -l app=lomad-power"
echo "  kubectl logs -n ${NAMESPACE} job/<job-name>"
