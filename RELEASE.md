# Release Guide for fragments-h5

This guide explains how to build and push Docker images and packages after making changes.

## Prerequisites

1. **Docker** installed and running
2. **GitHub CLI** (`gh`) installed and authenticated:
   ```bash
   gh auth login
   ```
3. **Conda** installed (for conda package builds)

## Current Version

The version is automatically read from `pyproject.toml` (currently **2.5.0**).

## Building and Pushing Docker Image

The Docker image will be pushed to `ghcr.io/nboley/fragments-h5:2.5.0` and `ghcr.io/nboley/fragments-h5:latest`.

### Quick Command
```bash
cd /home/nathanboley/src/fragments_h5
make push
```

This will:
1. Build the Docker image locally
2. Authenticate with GHCR using GitHub CLI
3. Tag the image with version and `latest`
4. Push both tags to GHCR

### Step-by-Step

1. **Build Docker image:**
   ```bash
   make docker
   ```
   This creates: `fragments-h5:2.5.0` and `fragments-h5:latest`

2. **Push to GHCR:**
   ```bash
   make push
   ```
   This authenticates, tags, and pushes to `ghcr.io/nboley/fragments-h5:2.5.0`

### Custom Configuration

Override defaults with environment variables:
```bash
GITHUB_USER=your-org make push  # Use different GitHub org/user
VERSION=2.5.0 make push         # Override version (defaults to pyproject.toml)
```

## Building Conda Package

### Build Only
```bash
make conda
```
Output: `conda-build-output/`

### Build and Upload to JFrog (if configured)
```bash
./scripts/build_conda_package.sh --upload
```

Requires environment variables:
- `JFROG_URL` - Artifactory URL
- `JFROG_REPO` - Repository name
- `JFROG_ACCESS_TOKEN` or `JFROG_USER`/`JFROG_PASSWORD`

## Creating Git Tag

After building and pushing, create a git tag:
```bash
make tag
```

This creates and pushes tag `v2.5.0` to the repository.

## Complete Release Workflow

To do everything at once:
```bash
# 1. Build conda package
make conda

# 2. Build and push Docker image
make push

# 3. Create and push git tag
make tag
```

Or use the `all` target (builds conda + docker + push, but not tag):
```bash
make all
make tag  # Still need to tag separately
```

## Verification

After pushing, verify the Docker image:
```bash
docker pull ghcr.io/nboley/fragments-h5:2.5.0
docker run --rm ghcr.io/nboley/fragments-h5:2.5.0 build-fragments-h5 --help
```

## Troubleshooting

- **Docker push fails**: Ensure `gh auth login` is completed
- **Version mismatch**: Check `pyproject.toml` and `conda-recipe/meta.yaml` both have version 2.5.0
- **Conda build fails**: Ensure conda-forge and bioconda channels are available
