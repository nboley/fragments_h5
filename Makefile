# fragments-h5 Makefile
#
# Usage:
#   make conda-build    # Build conda package
#   make conda          # Build and upload conda package
#   make docker-build   # Build Docker image
#   make docker-push    # Push Docker image to GHCR
#   make docker         # Build and push Docker image
#   make tag            # Create and push git tag
#   make all            # Build/upload conda, tag repo, build/push docker
#   make clean          # Remove build artifacts
#
# Configuration (override with environment variables):
#   GITHUB_USER    - GitHub username/org (default: nboley)
#   VERSION        - Package version (default: read from pyproject.toml)
#   DOCKER_PLATFORM - Platform for docker build (default: linux/amd64)

GITHUB_USER ?= nboley
VERSION ?= $(shell grep 'version = ' pyproject.toml | head -1 | sed 's/.*"\(.*\)".*/\1/')
IMAGE_NAME = fragments-h5
GHCR_IMAGE = ghcr.io/$(GITHUB_USER)/$(IMAGE_NAME)

.PHONY: all login conda-login docker-login require-clean-tree check-pyproject-clean conda-build conda-publish conda docker-build docker-push docker tag clean help

help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  login         Verify credentials for conda and docker"
	@echo "  conda-login   Verify JFrog credentials for conda publishing"
	@echo "  docker-login  Verify GitHub authentication for GHCR"
	@echo "  conda-build   Build conda package with rattler-build"
	@echo "  conda         Build and publish conda package"
	@echo "  docker-build  Build Docker image"
	@echo "  docker-push   Push Docker image to GHCR"
	@echo "  docker        Build and push Docker image"
	@echo "  tag           Create and push git tag v\$$VERSION"
	@echo "  all           Build/upload conda, tag repo, build/push docker"
	@echo "  clean         Remove build artifacts"
	@echo ""
	@echo "Configuration:"
	@echo "  GITHUB_USER=$(GITHUB_USER)"
	@echo "  VERSION=$(VERSION)"
	@echo "  GHCR_IMAGE=$(GHCR_IMAGE)"

all: login tag conda docker clean
	@echo ""
	@echo "========================================"
	@echo "Release $(VERSION) complete!"
	@echo "  ✓ Git tagged: v$(VERSION)"
	@echo "  ✓ Conda package built and published"
	@echo "  ✓ Docker pushed: $(GHCR_IMAGE):$(VERSION)"
	@echo "  ✓ Build artifacts cleaned"
	@echo "========================================"

docker-login:
	@echo "Verifying GitHub authentication..."
	@if ! command -v gh >/dev/null 2>&1; then \
		echo "❌ Error: GitHub CLI (gh) is not installed"; \
		exit 1; \
	fi
	@if ! gh auth status >/dev/null 2>&1; then \
		echo "❌ Error: Not authenticated with GitHub CLI. Run: gh auth login"; \
		exit 1; \
	fi
	@echo "✓ GitHub CLI authenticated"

login: conda-login docker-login
	@echo ""
	@echo "✓ All credentials verified successfully!"
	@echo ""

require-clean-tree:
	@# Whole-tree cleanliness gate for anything that produces a distributable
	@# artifact (git tag, conda package, Docker image). Tracked changes
	@# (staged or unstaged) and untracked files all disqualify. This is what
	@# prevents an artifact from disagreeing with the commit it claims to be
	@# built from -- the root cause of the v2.10.1 image/tag mismatch.
	@if ! git diff --quiet HEAD; then \
		echo "Error: working tree has uncommitted changes; refusing to proceed."; \
		echo "  Commit or stash your changes first."; \
		exit 1; \
	fi
	@if [ -n "$$(git ls-files --others --exclude-standard)" ]; then \
		echo "Error: working tree has untracked files; refusing to proceed."; \
		echo "  Add them to .gitignore or remove them first."; \
		exit 1; \
	fi

conda-build: require-clean-tree
	@echo "Building conda package with rattler-build..."
	@rattler-build build \
		--recipe conda-recipe/recipe.yaml \
		--output-dir conda-build-output \
		--channel conda-forge \
		--channel bioconda \
		--variant-config conda-recipe/variant_config.yaml \
		--variant pkg_version=$(VERSION); \
	BUILD_EXIT=$$?; \
	if [ $$BUILD_EXIT -ne 0 ] && { [ ! -d conda-build-output ] || [ -z "$$(find conda-build-output -name '*.conda' 2>/dev/null)" ]; }; then \
		echo "❌ Error: Conda build failed (exit code $$BUILD_EXIT)"; \
		exit $$BUILD_EXIT; \
	elif [ $$BUILD_EXIT -ne 0 ]; then \
		echo "⚠️  Warning: Build succeeded but cleanup failed (exit code $$BUILD_EXIT) - this is a known rattler-build issue"; \
	fi
	@echo "Conda package built: conda-build-output/"

conda: conda-build
	@echo "Conda package built successfully!"

docker-build: require-clean-tree
	@echo "Building Docker image $(IMAGE_NAME):$(VERSION)..."
	docker build \
		--build-arg BUILD_CODE_REVISION="$$(git describe --tags --always --dirty)" \
		-t $(IMAGE_NAME):$(VERSION) -t $(IMAGE_NAME):latest .
	@echo "Docker image built: $(IMAGE_NAME):$(VERSION)"

docker-push: docker-build
	@echo "Authenticating with GHCR..."
	@docker logout ghcr.io 2>/dev/null || true
	@gh auth token | docker login ghcr.io -u $(GITHUB_USER) --password-stdin
	@echo "Tagging images..."
	docker tag $(IMAGE_NAME):$(VERSION) $(GHCR_IMAGE):$(VERSION)
	docker tag $(IMAGE_NAME):latest $(GHCR_IMAGE):latest
	@echo "Pushing to GHCR..."
	docker push $(GHCR_IMAGE):$(VERSION)
	docker push $(GHCR_IMAGE):latest
	@echo "Pushed: $(GHCR_IMAGE):$(VERSION)"

docker: docker-build docker-push
	@echo "Docker image built and pushed successfully!"

check-pyproject-clean:
	@# VERSION is read from the WORKING TREE. If pyproject.toml is uncommitted,
	@# the tag would point at HEAD, which declares a different version. This is
	@# how v2.10.1 came to point at a commit whose pyproject.toml said 2.10.0.
	@# This check runs before the whole-tree require-clean-tree gate so that,
	@# in the common case (only pyproject.toml is dirty), the error explains
	@# *why* that specifically matters for tagging rather than just refusing.
	@if ! git diff --quiet HEAD -- pyproject.toml; then \
		echo "Error: pyproject.toml has uncommitted changes; refusing to tag."; \
		echo "  Working tree declares $(VERSION), but the tag would point at HEAD,"; \
		echo "  which declares something else. Commit the version bump first."; \
		exit 1; \
	fi

tag: check-pyproject-clean require-clean-tree
	@if git rev-parse "v$(VERSION)" >/dev/null 2>&1; then \
		echo "Error: Tag v$(VERSION) already exists"; \
		echo "  Bump the version in pyproject.toml first."; \
		exit 1; \
	fi
	@echo "Creating git tag v$(VERSION)..."
	git tag -a v$(VERSION) -m "Release v$(VERSION)"
	@echo "Pushing tag to origin..."
	git push origin v$(VERSION)
	@echo "Tagged: v$(VERSION)"

clean:
	rm -rf conda-build-output/
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf src/*.egg-info/
	rm -f src/fragments_h5/sequence.c
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	@echo "Cleaned build artifacts"

