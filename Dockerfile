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
ARG BUILD_CODE_REVISION=""
RUN if [ -n "$BUILD_CODE_REVISION" ]; then \
        printf 'REVISION = "%s"\n' "$BUILD_CODE_REVISION" \
            > src/fragments_h5/_build_revision.py; \
    fi && \
    python setup.py build_ext --inplace && \
    pip install --no-cache-dir --no-deps . && \
    # awscrt is installed explicitly because the line above uses --no-deps, so
    # nothing in pyproject.toml reaches this image. repair.py uploads with
    # ChecksumAlgorithm="CRC64NVME", which botocore can only compute with the
    # AWS CRT bindings; the conda env does not provide them. Installed with pip
    # rather than conda because the cp313-abi3 wheel works on this image's
    # Python 3.14, whereas the conda-forge builds lag the interpreter.
    pip install --no-cache-dir awscrt && \
    # Fail the BUILD rather than the 218th production upload. This assertion
    # exists because the missing dependency was invisible until the write path
    # ran for the first time against real S3, three minutes into a repair.
    python -c "import awscrt, botocore; print('awscrt', awscrt.__version__)"

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
# Install s5cmd for high-performance S3 staging with Nextflow
# Install aws-cli v2 (required by Nextflow for bin folder sync)
ARG S5CMD_VERSION=2.3.0
RUN apt-get update && apt-get install -y --no-install-recommends procps curl unzip && \
    curl -sL https://github.com/peak/s5cmd/releases/download/v${S5CMD_VERSION}/s5cmd_${S5CMD_VERSION}_Linux-64bit.tar.gz | tar xz -C /usr/local/bin/ && \
    chmod +x /usr/local/bin/s5cmd && \
    curl -sL "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "/tmp/awscliv2.zip" && \
    unzip -q /tmp/awscliv2.zip -d /tmp && \
    /tmp/aws/install && \
    rm -rf /tmp/awscliv2.zip /tmp/aws && \
    apt-get remove -y curl unzip && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

# Copy only the cleaned conda environment
COPY --from=builder /opt/conda /opt/conda

WORKDIR /data
ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["build-fragments-h5", "--help"]
