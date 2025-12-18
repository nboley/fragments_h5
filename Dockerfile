# Using micromamba for fast, small conda-based builds
FROM mambaorg/micromamba:1.5-bookworm-slim AS builder

USER root

# Create environment with dependencies
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Copy and install the package
COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp/fragments_h5
WORKDIR /tmp/fragments_h5

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN python setup.py build_ext --inplace && \
    pip install --no-cache-dir --no-deps .

# Remove build-time-only packages (but not pip which would remove python)
RUN micromamba remove -y -n base c-compiler cython && \
    micromamba clean --all --yes && \
    find /opt/conda -name "*.a" -delete && \
    find /opt/conda -name "*.pyc" -delete && \
    find /opt/conda -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true && \
    find /opt/conda -name "tests" -type d -exec rm -rf {} + 2>/dev/null || true && \
    find /opt/conda -name "test" -type d -exec rm -rf {} + 2>/dev/null || true && \
    find /opt/conda -name "*.pyx" -delete && \
    rm -rf /opt/conda/pkgs/* && \
    rm -rf /opt/conda/share/man && \
    rm -rf /opt/conda/share/doc && \
    rm -rf /opt/conda/share/gtk-doc

# Runtime stage
FROM mambaorg/micromamba:1.5-bookworm-slim

USER root

# Install procps (provides 'ps' command needed by Nextflow for metrics)
RUN apt-get update && apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*

# Copy AWS CLI v2 from official Amazon image
COPY --from=amazon/aws-cli:latest /usr/local/aws-cli /usr/local/aws-cli
RUN ln -s /usr/local/aws-cli/v2/current/bin/aws /usr/local/bin/aws

# Copy only the cleaned conda environment
COPY --from=builder /opt/conda /opt/conda

WORKDIR /data
ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["build-fragments-h5", "--help"]
