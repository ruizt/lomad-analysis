#!/bin/bash
# build-image.sh — build and push the shared simulation image
#
# Usage (from the repo root):
#   bash simulations/build-image.sh              # lomad from GitHub, push
#   bash simulations/build-image.sh --local      # lomad from ../lomad-package
#   bash simulations/build-image.sh --no-push    # build only, do not push
#
# By default lomad is installed from GitHub at LOMAD_REF (default: main),
# which requires ruizt/lomad-package to be public.
#
# --local instead builds lomad from the sibling checkout and installs that
# source tree. Use it while the package repo is private, or to build an image
# matching your working tree rather than a branch -- which is what you want
# when re-running simulations against an unreleased change.
#
# Requires Docker running, and `docker login ghcr.io` for pushing.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PKG_DIR="${LOMAD_PKG_DIR:-${REPO_ROOT}/../lomad-package}"
BUILD_DIR="${REPO_ROOT}/simulations/_build"
IMAGE="${IMAGE:-ghcr.io/ruizt/lomad-simulations:latest}"

PUSH="--push"
SOURCE="remote"
LOMAD_REF="${LOMAD_REF:-main}"
for arg in "$@"; do
  case "$arg" in
    --no-push) PUSH="--load" ;;
    --local)   SOURCE="local" ;;
    *) echo "unknown option: $arg" >&2; exit 1 ;;
  esac
done

# ---- Prepare the lomad source ------------------------------------------------
# _build must exist either way: the Dockerfile COPYs it unconditionally and
# only reads it when LOMAD_SOURCE=local.

rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"

if [ "${SOURCE}" = "local" ]; then
  [ -f "${PKG_DIR}/DESCRIPTION" ] || {
    echo "error: no R package at ${PKG_DIR}" >&2
    echo "       expected lomad-package as a sibling of this repo," >&2
    echo "       or set LOMAD_PKG_DIR to its location." >&2
    exit 1
  }

  PKG_VERSION=$(awk '/^Version:/ {print $2}' "${PKG_DIR}/DESCRIPTION")
  PKG_SHA=$(git -C "${PKG_DIR}" rev-parse --short HEAD 2>/dev/null || echo "unknown")
  echo "lomad ${PKG_VERSION} (${PKG_SHA}) from ${PKG_DIR}"

  if [ -n "$(git -C "${PKG_DIR}" status --porcelain 2>/dev/null | grep -v '^??' || true)" ]; then
    echo "note: lomad-package has uncommitted changes; they will go into the image."
  fi

  # R CMD build first so .Rbuildignore is honoured, then unpack: pak's local::
  # ref resolves a package directory, not an archive.
  ( cd "${BUILD_DIR}" && R CMD build --no-build-vignettes --no-manual "${PKG_DIR}" >/dev/null )
  tar -xzf "${BUILD_DIR}"/lomad_*.tar.gz -C "${BUILD_DIR}"
  rm -f "${BUILD_DIR}"/lomad_*.tar.gz
  echo "Source tree: $(du -sh "${BUILD_DIR}/lomad" | cut -f1)"
else
  PKG_VERSION="remote"
  PKG_SHA="${LOMAD_REF}"
  echo "lomad from ruizt/lomad-package@${LOMAD_REF} (requires the repo to be public)"
  touch "${BUILD_DIR}/.unused"
fi

# ---- Build the image ---------------------------------------------------------

echo "Building ${IMAGE} ..."
docker buildx build --platform linux/amd64 \
  -f "${REPO_ROOT}/simulations/Dockerfile" \
  -t "${IMAGE}" \
  --build-arg "LOMAD_SOURCE=${SOURCE}" \
  --build-arg "LOMAD_REF=${LOMAD_REF}" \
  --label "org.opencontainers.image.source=https://github.com/ruizt/lomad-analysis" \
  --label "lomad.version=${PKG_VERSION}" \
  --label "lomad.revision=${PKG_SHA}" \
  ${PUSH} "${REPO_ROOT}"

echo ""
echo "Done: ${IMAGE}  (lomad ${PKG_VERSION}, ${PKG_SHA}, source=${SOURCE})"
if [ "${PUSH}" = "--push" ]; then
  echo "The ghcr package must be Public for the cluster to pull it."
fi
