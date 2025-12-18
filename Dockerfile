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

# Remove build-time-only packages and cleanup to reduce image size
# Note: removing pip/setuptools via micromamba cascades to remove python, so we remove manually
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
    rm -rf /opt/conda/share/gtk-doc && \
    rm -rf /opt/conda/share/terminfo && \
    rm -rf /opt/conda/share/locale && \
    rm -rf /opt/conda/lib/python*/ensurepip && \
    rm -rf /opt/conda/lib/python*/idlelib && \
    rm -rf /opt/conda/lib/python*/tkinter && \
    rm -rf /opt/conda/lib/python*/lib2to3 && \
    rm -rf /opt/conda/lib/python*/site-packages/pip* && \
    rm -rf /opt/conda/lib/python*/site-packages/setuptools* && \
    rm -rf /opt/conda/lib/python*/site-packages/pkg_resources && \
    rm -rf /opt/conda/lib/python*/site-packages/_distutils_hack && \
    rm -rf /opt/conda/lib/python*/site-packages/distutils-precedence.pth && \
    rm -rf /opt/conda/bin/pip* && \
    rm -rf /opt/conda/include

# Runtime stage
FROM mambaorg/micromamba:1.5-bookworm-slim

USER root

# Install procps (provides 'ps' command needed by Nextflow for metrics)
RUN apt-get update && apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*

# Install s5cmd (fast S3 client, ~20MB vs ~250MB for aws-cli)
ARG S5CMD_VERSION=2.2.2
ARG TARGETARCH
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates && \
    ARCH=$([ "$TARGETARCH" = "arm64" ] && echo "arm64" || echo "amd64") && \
    curl -sL "https://github.com/peak/s5cmd/releases/download/v${S5CMD_VERSION}/s5cmd_${S5CMD_VERSION}_Linux-${ARCH}.tar.gz" | \
    tar xz -C /usr/local/bin s5cmd && \
    chmod +x /usr/local/bin/s5cmd && \
    apt-get remove -y curl && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

# Create 'aws' shim that translates 'aws s3' commands to s5cmd
RUN printf '%s\n' '#!/bin/bash' \
    'if [[ "$1" == "s3" ]]; then' \
    '    shift' \
    '    exec /usr/local/bin/s5cmd "$@"' \
    'fi' \
    'echo "Unsupported aws command: $*" >&2' \
    'exit 1' > /usr/local/bin/aws && chmod +x /usr/local/bin/aws

# Copy only the cleaned conda environment
COPY --from=builder /opt/conda /opt/conda

WORKDIR /data
ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["build-fragments-h5", "--help"]
