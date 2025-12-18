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

# Remove build-time-only packages and aggressive cleanup to reduce image size
RUN micromamba remove -y -n base c-compiler cython setuptools pip && \
    micromamba clean --all --yes && \
    find /opt/conda -name "*.a" -delete && \
    find /opt/conda -name "*.pyc" -delete && \
    find /opt/conda -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true && \
    find /opt/conda -name "tests" -type d -exec rm -rf {} + 2>/dev/null || true && \
    find /opt/conda -name "test" -type d -exec rm -rf {} + 2>/dev/null || true && \
    find /opt/conda -name "*.pyx" -delete && \
    find /opt/conda -name "*.c" -delete && \
    find /opt/conda -name "*.h" ! -path "*/numpy/*" -delete 2>/dev/null || true && \
    rm -rf /opt/conda/pkgs/* && \
    rm -rf /opt/conda/share/man && \
    rm -rf /opt/conda/share/doc && \
    rm -rf /opt/conda/share/gtk-doc && \
    rm -rf /opt/conda/include && \
    rm -rf /opt/conda/conda-meta/*.json && \
    strip --strip-unneeded /opt/conda/lib/*.so* 2>/dev/null || true

# Runtime stage
FROM mambaorg/micromamba:1.5-bookworm-slim

# Copy only the cleaned conda environment
COPY --from=builder /opt/conda /opt/conda

WORKDIR /data
ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["build-fragments-h5", "--help"]
