#!/bin/bash
# submit_sweep.sh — submit one job per value of d
# Usage: bash lomad-sim/submit_sweep.sh

D_VALUES=(0.5 1 2 3 5)

for d in "${D_VALUES[@]}"; do
  # Create a job name safe for Kubernetes (replace . with -)
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
  template:
    spec:
      restartPolicy: Never
      imagePullSecrets:
        - name: ghcr-secret
      containers:
        - name: lomad-sim
          image: ghcr.io/otishunt/lomad-sim:latest
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
              value: "30"
            - name: SIM_Q
              value: "30"
            - name: SIM_B
              value: "1000"
            - name: SIM_METHOD
              value: "boot"
            - name: SIM_NCORES
              value: "1"
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

echo "All jobs submitted. Monitor with: kubectl get jobs"
