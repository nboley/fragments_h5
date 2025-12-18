# fragments-h5 Makefile
#
# Usage:
#   make conda          # Build conda package
#   make docker         # Build Docker image
#   make push           # Push Docker image to GHCR
#   make all            # Build conda + docker + push
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

.PHONY: all conda docker push tag clean help

help:
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  conda    Build conda package"
	@echo "  docker   Build Docker image"
	@echo "  push     Tag and push Docker image to GHCR"
	@echo "  tag      Create and push git tag v\$$VERSION"
	@echo "  all      Build everything and push"
	@echo "  clean    Remove build artifacts"
	@echo ""
	@echo "Configuration:"
	@echo "  GITHUB_USER=$(GITHUB_USER)"
	@echo "  VERSION=$(VERSION)"
	@echo "  GHCR_IMAGE=$(GHCR_IMAGE)"

all: conda docker push

conda:
	@echo "Building conda package..."
	conda build conda-recipe \
		--output-folder conda-build-output \
		--channel conda-forge \
		--channel bioconda
	@echo "Conda package built: conda-build-output/"

docker:
	@echo "Building Docker image $(IMAGE_NAME):$(VERSION)..."
	docker build -t $(IMAGE_NAME):$(VERSION) -t $(IMAGE_NAME):latest .
	@echo "Docker image built: $(IMAGE_NAME):$(VERSION)"

push: docker
	@echo "Authenticating with GHCR..."
	@docker logout ghcr.io 2>/dev/null || true
	@gh auth token | docker login ghcr.io -u $(GITHUB_USER) --password-stdin
	@echo "Tagging for GHCR..."
	docker tag $(IMAGE_NAME):$(VERSION) $(GHCR_IMAGE):$(VERSION)
	docker tag $(IMAGE_NAME):latest $(GHCR_IMAGE):latest
	@echo "Pushing to GHCR..."
	docker push $(GHCR_IMAGE):$(VERSION)
	docker push $(GHCR_IMAGE):latest
	@echo "Pushed: $(GHCR_IMAGE):$(VERSION)"

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
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	@echo "Cleaned build artifacts"

