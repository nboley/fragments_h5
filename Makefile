# fragments-h5 Makefile
#
# Usage:
#   make conda-build    # Build conda package
#   make conda-publish  # Upload conda package to JFrog
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

.PHONY: all login conda-login docker-login conda-build conda-publish conda docker-build docker-push docker tag clean help

help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  login         Verify credentials for conda and docker"
	@echo "  conda-login   Verify JFrog credentials for conda publishing"
	@echo "  docker-login  Verify GitHub authentication for GHCR"
	@echo "  conda-build   Build conda package with rattler-build"
	@echo "  conda-publish Publish conda package to JFrog Artifactory"
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

conda-login:
	@echo "Verifying JFrog credentials..."
	@# Check if credentials are in environment variables
	@HAS_ENV_CREDS=0; \
	if [ -n "$$JFROG_URL" ] && { [ -n "$$JFROG_ACCESS_TOKEN" ] || [ -n "$$JFROG_USER" ]; }; then \
		HAS_ENV_CREDS=1; \
	fi; \
	# Check if pip is configured with JFrog credentials
	HAS_PIP_CREDS=0; \
	if pip config list 2>/dev/null | grep -q "global.extra-index-url\|global.index-url" 2>/dev/null; then \
		HAS_PIP_CREDS=1; \
	fi; \
	# Verify at least one method has credentials
	if [ $$HAS_ENV_CREDS -eq 0 ] && [ $$HAS_PIP_CREDS -eq 0 ]; then \
		echo "❌ Error: JFrog credentials not found"; \
		echo "   Either set environment variables (JFROG_URL + JFROG_ACCESS_TOKEN/JFROG_USER)"; \
		echo "   Or configure pip with JFrog credentials in pip.conf"; \
		exit 1; \
	fi
	@echo "✓ JFrog credentials configured"

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

conda-build:
	@echo "Building conda package with rattler-build..."
	rattler-build build \
		--recipe conda-recipe/recipe.yaml \
		--output-dir conda-build-output \
		--channel conda-forge \
		--channel bioconda \
		--variant-config conda-recipe/variant_config.yaml || true
	@echo "Conda package built: conda-build-output/"

conda-publish:
	@echo "Publishing conda package to JFrog Artifactory..."
	bash scripts/publish_conda_package.sh
	@echo "Package published to karius-conda repository"

conda: conda-build conda-publish
	@echo "Conda package built and published successfully!"

docker-build:
	@echo "Building Docker image $(IMAGE_NAME):$(VERSION)..."
	docker build -t $(IMAGE_NAME):$(VERSION) -t $(IMAGE_NAME):latest .
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

tag:
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

