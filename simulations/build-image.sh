#!/bin/bash
# build-image.sh — build and push the shared simulation image
#
# Usage (from the repo root):
#   bash simulations/build-image.sh            # build and push :latest
#   bash simulations/build-image.sh --no-push  # build locally only
#
# Builds the lomad source tarball from the sibling lomad-package checkout and
# installs it into the image from the build context. This is deliberate: the
# package repo is private, so pak inside the container cannot authenticate to
# GitHub. Building locally also means the image always matches the working
# tree, which is what you want when re-running simulations against a change
# that has not been released yet.
#
# Requires Docker running, and `docker login ghcr.io` for pushing.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PKG_DIR="${LOMAD_PKG_DIR:-${REPO_ROOT}/../lomad-package}"
BUILD_DIR="${REPO_ROOT}/simulations/_build"
IMAGE="${IMAGE:-ghcr.io/ruizt/lomad-sims:latest}"

PUSH="--push"
[ "${1:-}" = "--no-push" ] && PUSH="--load"

# ---- Locate the package ------------------------------------------------------

[ -f "${PKG_DIR}/DESCRIPTION" ] || {
  echo "error: no R package at ${PKG_DIR}" >&2
  echo "       expected lomad-package as a sibling of this repo," >&2
  echo "       or set LOMAD_PKG_DIR to its location." >&2
  exit 1
}

PKG_VERSION=$(awk '/^Version:/ {print $2}' "${PKG_DIR}/DESCRIPTION")
PKG_SHA=$(git -C "${PKG_DIR}" rev-parse --short HEAD 2>/dev/null || echo "unknown")
echo "Building lomad ${PKG_VERSION} (${PKG_SHA}) from ${PKG_DIR}"

if [ -n "$(git -C "${PKG_DIR}" status --porcelain 2>/dev/null | grep -v '^??' || true)" ]; then
  echo "note: lomad-package has uncommitted changes; they will go into the image."
fi

# ---- Build the source tarball ------------------------------------------------

rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"

# R CMD build first so .Rbuildignore is honoured, then unpack: pak's local::
# ref resolves a package directory, not an archive.
( cd "${BUILD_DIR}" && R CMD build --no-build-vignettes --no-manual "${PKG_DIR}" >/dev/null )
tar -xzf "${BUILD_DIR}"/lomad_*.tar.gz -C "${BUILD_DIR}"
rm -f "${BUILD_DIR}"/lomad_*.tar.gz
echo "Source tree: $(du -sh "${BUILD_DIR}/lomad" | cut -f1)"

# ---- Build the image ---------------------------------------------------------

echo "Building ${IMAGE} ..."
docker buildx build --platform linux/amd64 \
  -f "${REPO_ROOT}/simulations/Dockerfile" \
  -t "${IMAGE}" \
  --label "org.opencontainers.image.source=https://github.com/ruizt/lomad-analysis" \
  --label "lomad.version=${PKG_VERSION}" \
  --label "lomad.revision=${PKG_SHA}" \
  ${PUSH} "${REPO_ROOT}"

echo ""
echo "Done: ${IMAGE}  (lomad ${PKG_VERSION}, ${PKG_SHA})"
if [ "${PUSH}" = "--push" ]; then
  echo "The ghcr package must be Public for the cluster to pull it."
fi
