#!/bin/bash
#
# Build and upload fragments-h5 conda package to JFrog Artifactory
#
# Usage:
#   ./scripts/build_conda_package.sh                    # Build only
#   ./scripts/build_conda_package.sh --upload           # Build and upload
#
# Environment variables:
#   JFROG_URL      - JFrog Artifactory URL (e.g., https://company.jfrog.io/artifactory)
#   JFROG_REPO     - Repository name (e.g., conda-local)
#   JFROG_USER     - JFrog username (or use JFROG_ACCESS_TOKEN)
#   JFROG_PASSWORD - JFrog password (or use JFROG_ACCESS_TOKEN)
#   JFROG_ACCESS_TOKEN - JFrog access token (alternative to user/password)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

UPLOAD=false
[[ "${1:-}" == "--upload" ]] && UPLOAD=true

# Build
echo "Building conda package..."
conda build "$PROJECT_ROOT/conda-recipe" \
    --output-folder "$PROJECT_ROOT/conda-build-output" \
    --channel conda-forge \
    --channel bioconda

PACKAGE=$(conda build "$PROJECT_ROOT/conda-recipe" --output 2>/dev/null | tail -1)
echo "Built: $PACKAGE"

# Upload to JFrog
if [ "$UPLOAD" = true ]; then
    : "${JFROG_URL:?Set JFROG_URL environment variable}"
    : "${JFROG_REPO:?Set JFROG_REPO environment variable}"
    
    # Determine auth method
    if [ -n "${JFROG_ACCESS_TOKEN:-}" ]; then
        AUTH="-H \"Authorization: Bearer $JFROG_ACCESS_TOKEN\""
    else
        : "${JFROG_USER:?Set JFROG_USER or JFROG_ACCESS_TOKEN}"
        : "${JFROG_PASSWORD:?Set JFROG_PASSWORD or JFROG_ACCESS_TOKEN}"
        AUTH="-u $JFROG_USER:$JFROG_PASSWORD"
    fi
    
    FILENAME=$(basename "$PACKAGE")
    SUBDIR=$(basename "$(dirname "$PACKAGE")")  # e.g., "noarch" or "linux-64"
    
    echo "Uploading to JFrog..."
    eval curl -f $AUTH -T "$PACKAGE" \
        "$JFROG_URL/$JFROG_REPO/$SUBDIR/$FILENAME"
    
    echo "Done! Install with: conda install -c $JFROG_URL/$JFROG_REPO fragments-h5"
fi
