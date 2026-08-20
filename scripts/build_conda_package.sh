#!/bin/bash
#
# Build fragments-h5 conda package using rattler-build
#
# Usage:
#   ./scripts/build_conda_package.sh
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Build
echo "Building conda package with rattler-build..."
VERSION=$(grep '^version' "$PROJECT_ROOT/pyproject.toml" | head -1 | sed -E 's/.*"([0-9.]+)".*/\1/')
rattler-build build \
    --recipe "$PROJECT_ROOT/conda-recipe/recipe.yaml" \
    --output-dir "$PROJECT_ROOT/conda-build-output" \
    --channel conda-forge \
    --channel bioconda \
    --variant-config "$PROJECT_ROOT/conda-recipe/variant_config.yaml" \
    --variant "pkg_version=${VERSION}"

# Get the built package path
PACKAGE=$(find "$PROJECT_ROOT/conda-build-output" -name "fragments-h5-*.conda" -o -name "fragments-h5-*.tar.bz2" | head -1)
echo "Built: $PACKAGE"
echo ""
echo "To publish to JFrog Artifactory, run:"
echo "  make conda-publish"
echo "Or to build and publish in one step:"
echo "  make conda"

