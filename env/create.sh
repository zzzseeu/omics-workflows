#!/usr/bin/env bash
# Create once; do not update an existing public environment in place.
set -Eeuo pipefail
REPO_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
MAMBA=${MICROMAMBA_EXE:-micromamba}
ENV_PREFIX=${ASSEMBLY_ENV_PREFIX:-${OMICS_ENV_PREFIX:-}}
: "${ENV_PREFIX:?Set OMICS_ENV_PREFIX or ASSEMBLY_ENV_PREFIX to the unified environment prefix}"
RESOURCE_ROOT=${OMICS_RESOURCE_ROOT:-}
DEPLOYMENT_DIR=${ASSEMBLY_DEPLOYMENT_DIR:-}
if [[ -z "$DEPLOYMENT_DIR" ]]; then
    : "${RESOURCE_ROOT:?Set ASSEMBLY_DEPLOYMENT_DIR or OMICS_RESOURCE_ROOT outside the code repository}"
    DEPLOYMENT_DIR="$RESOURCE_ROOT/deployments/assembly_$(date -u +%Y%m%dT%H%M%SZ)"
fi
DEFAULT_SPEC="$REPO_ROOT/env/environment.yaml"
# linux-64.lock records the historical environment; use the complete seven-stage
# dependency specification unless a newly resolved explicit lock is supplied.
ENV_SPEC=${ASSEMBLY_ENV_SPEC:-$DEFAULT_SPEC}
mkdir -p "$DEPLOYMENT_DIR"
if [[ -e "$ENV_PREFIX" ]]; then
    echo "Refusing to overwrite existing environment: $ENV_PREFIX" >&2
    exit 2
fi
trap 'rc=$?; if ((rc != 0)); then printf "FAIL\texit=%s\n" "$rc" > "$DEPLOYMENT_DIR/status.tsv"; fi' EXIT
printf 'RUNNING\t%s\n' "$(date -u --iso-8601=seconds)" > "$DEPLOYMENT_DIR/status.tsv"
cp "$REPO_ROOT/env/environment.yaml" "$DEPLOYMENT_DIR/assembly.requested.yaml"
cp "$ENV_SPEC" "$DEPLOYMENT_DIR/installation.spec"
"$MAMBA" --version > "$DEPLOYMENT_DIR/micromamba.version.txt"
"$MAMBA" info --json > "$DEPLOYMENT_DIR/micromamba.info.json"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# Avoid thousands of small HTTP requests through the cluster's network proxy.
export MAMBA_USE_SHARDED_REPODATA=${MAMBA_USE_SHARDED_REPODATA:-false}
"$MAMBA" create --yes --prefix "$ENV_PREFIX" \
    --override-channels --strict-channel-priority \
    --channel https://conda.anaconda.org/conda-forge \
    --channel https://conda.anaconda.org/bioconda \
    --file "$ENV_SPEC" \
    > "$DEPLOYMENT_DIR/create.log" 2>&1
"$ENV_PREFIX/bin/python" "$REPO_ROOT/env/install_haphic.py" \
    --env-prefix "$ENV_PREFIX" > "$DEPLOYMENT_DIR/haphic.install.log" 2>&1
"$MAMBA" list --prefix "$ENV_PREFIX" --json > "$DEPLOYMENT_DIR/packages.json"
"$MAMBA" env export --prefix "$ENV_PREFIX" --explicit > "$DEPLOYMENT_DIR/assembly.linux-64.pin.txt"
"$MAMBA" env export --prefix "$ENV_PREFIX" > "$DEPLOYMENT_DIR/assembly.resolved.yaml"
if [[ ${ASSEMBLY_VALIDATE:-1} == 0 ]]; then
    printf 'INSTALLED\t%s\n' "$(date -u --iso-8601=seconds)" > "$DEPLOYMENT_DIR/status.tsv"
    printf 'prefix=%s\ndeployment=%s\n' "$ENV_PREFIX" "$DEPLOYMENT_DIR" > "$DEPLOYMENT_DIR/INSTALLED.ok"
    echo "Environment installed; validation is pending: $DEPLOYMENT_DIR"
    exit 0
fi
"$MAMBA" run --prefix "$ENV_PREFIX" \
    python "$REPO_ROOT/env/validate.py" \
    --output "$DEPLOYMENT_DIR/validation" > "$DEPLOYMENT_DIR/validation.log" 2>&1
sha256sum "$DEPLOYMENT_DIR/assembly.linux-64.pin.txt" > "$DEPLOYMENT_DIR/lock.sha256"
printf 'PASS\t%s\n' "$(date -u --iso-8601=seconds)" > "$DEPLOYMENT_DIR/status.tsv"
printf 'ASSEMBLY_ENVIRONMENT=PASS\nprefix=%s\ndeployment=%s\n' "$ENV_PREFIX" "$DEPLOYMENT_DIR" > "$DEPLOYMENT_DIR/COMPLETED.ok"
cat "$DEPLOYMENT_DIR/COMPLETED.ok"
