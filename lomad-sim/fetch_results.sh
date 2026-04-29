#!/bin/bash
# fetch_results.sh — copy the lomad results CSV off the cluster PVC
#
# The simulation jobs write results to a Kubernetes PersistentVolumeClaim (PVC)
# that only pods on the cluster can access directly. This script:
#   1. Spins up a temporary busybox pod that mounts the PVC
#   2. Copies lomad_results_log.csv into lomad-sim/results/ in the repo
#   3. Cleans up the temporary pod immediately
#   4. Commits and pushes the updated CSV to GitHub automatically
#
# Usage (run from anywhere inside the repo):
#   bash lomad-sim/fetch_results.sh
#
# Skip the auto-push if you just want the file locally:
#   NO_PUSH=1 bash lomad-sim/fetch_results.sh

set -euo pipefail

# ── CONFIG ───────────────────────────────────────────────────────────────────
NAMESPACE="${NAMESPACE:-cal-poly-lomad}"
PVC_NAME="${PVC_NAME:-lomad-sim-results}"
REMOTE_FILE="/output/lomad_results_log.csv"
# Always save to a fixed path in the repo so the viewer artifact can find it
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOCAL_OUT="${REPO_ROOT}/lomad-sim/results/lomad_results_log.csv"
HELPER_POD="lomad-fetch-$(date +%s)"   # unique name avoids conflicts
# ─────────────────────────────────────────────────────────────────────────────

echo "Namespace : ${NAMESPACE}"
echo "PVC       : ${PVC_NAME}"
echo "Remote    : ${REMOTE_FILE}"
echo "Local     : ${LOCAL_OUT}"
echo ""

# Step 1: Start the helper pod (mounts the PVC read-only, sleeps long enough
# for us to copy the file, then the pod exits on its own after 60 s).
# We use 'sleep 60' as a minimal idle — the pod exits after the copy, so it
# never holds resources beyond the ~5 seconds it takes to copy the file.
echo "Starting helper pod ${HELPER_POD}..."
kubectl run "${HELPER_POD}" \
  --namespace="${NAMESPACE}" \
  --image=busybox \
  --restart=Never \
  --overrides="{
    \"spec\": {
      \"volumes\": [{
        \"name\": \"output\",
        \"persistentVolumeClaim\": {\"claimName\": \"${PVC_NAME}\"}
      }],
      \"containers\": [{
        \"name\": \"fetch\",
        \"image\": \"busybox\",
        \"command\": [\"sleep\", \"60\"],
        \"volumeMounts\": [{
          \"name\": \"output\",
          \"mountPath\": \"/output\",
          \"readOnly\": true
        }]
      }]
    }
  }"

# Step 2: Wait until the pod is Running
echo "Waiting for pod to be ready..."
kubectl wait pod "${HELPER_POD}" \
  --namespace="${NAMESPACE}" \
  --for=condition=Ready \
  --timeout=60s

# Step 3: Copy the CSV
echo "Copying results..."
kubectl cp \
  --namespace="${NAMESPACE}" \
  "${HELPER_POD}:${REMOTE_FILE}" \
  "${LOCAL_OUT}"

echo "Saved -> ${LOCAL_OUT}"

# Step 4: Delete the helper pod immediately — don't hold the node
echo "Cleaning up helper pod..."
kubectl delete pod "${HELPER_POD}" \
  --namespace="${NAMESPACE}" \
  --grace-period=0 \
  --ignore-not-found

echo ""

# Step 5: Commit and push to GitHub so results are visible without a terminal.
# Skip this step by setting NO_PUSH=1.
if [[ "${NO_PUSH:-0}" == "1" ]]; then
  echo "Skipping git push (NO_PUSH=1)."
  echo "File is at: ${LOCAL_OUT}"
else
  echo "Committing and pushing results to GitHub..."
  cd "${REPO_ROOT}"
  git add "${LOCAL_OUT}"
  # Only commit if there are actual changes (avoids empty commits)
  if git diff --cached --quiet; then
    echo "No new results to commit (CSV unchanged)."
  else
    ROW_COUNT=$(( $(wc -l < "${LOCAL_OUT}") - 1 ))
    git commit -m "results: update lomad_results_log.csv (${ROW_COUNT} rows) [$(date '+%Y-%m-%d %H:%M')]"
    git push
    echo "Pushed -> github.com/otishunt/lomad"
  fi
fi

echo ""
echo "Done. Open the lomad Results Viewer in Cowork to see the table."
