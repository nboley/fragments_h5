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
    rm -rf /opt/conda/include && \
    rm -rf /opt/conda/bin/*-ld && \
    rm -rf /opt/conda/bin/nghttp* && \
    rm -rf /opt/conda/bin/h5* && \
    rm -rf /opt/conda/bin/bzip2* && \
    rm -rf /opt/conda/bin/bunzip2 && \
    rm -rf /opt/conda/bin/bzcat && \
    rm -rf /opt/conda/bin/lz4* && \
    rm -rf /opt/conda/bin/zstd* && \
    rm -rf /opt/conda/man && \
    rm -rf /opt/conda/sbin && \
    rm -rf /opt/conda/conda-meta && \
    rm -rf /opt/conda/lib/tcl* && \
    rm -rf /opt/conda/lib/tk* && \
    rm -rf /opt/conda/lib/libtcl* && \
    rm -rf /opt/conda/lib/libtk* && \
    rm -rf /opt/conda/lib/libsqlite* && \
    rm -rf /opt/conda/lib/itcl* && \
    rm -rf /opt/conda/lib/tdbc* && \
    rm -rf /opt/conda/lib/thread* && \
    rm -rf /opt/conda/lib/sqlite* && \
    rm -rf /opt/conda/lib/libhdf5_cpp* && \
    rm -rf /opt/conda/lib/libhdf5_fortran* && \
    rm -rf /opt/conda/lib/libhdf5_hl_cpp* && \
    rm -rf /opt/conda/lib/libhdf5_hl_fortran* && \
    rm -rf /opt/conda/lib/cmake && \
    rm -rf /opt/conda/lib/pkgconfig && \
    rm -rf /opt/conda/lib/krb5/plugins && \
    strip --strip-unneeded /opt/conda/lib/*.so* 2>/dev/null || true && \
    strip --strip-unneeded /opt/conda/lib/python*/lib-dynload/*.so* 2>/dev/null || true

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
    ARCH=$([ "$TARGETARCH" = "arm64" ] && echo "arm64" || echo "64bit") && \
    curl -sL "https://github.com/peak/s5cmd/releases/download/v${S5CMD_VERSION}/s5cmd_${S5CMD_VERSION}_Linux-${ARCH}.tar.gz" | \
    tar xz -C /usr/local/bin s5cmd && \
    chmod +x /usr/local/bin/s5cmd && \
    apt-get remove -y curl && apt-get autoremove -y && \
    rm -rf /var/lib/apt/lists/*

# Create 'aws' shim that translates 'aws s3' commands to s5cmd
# Handles: aws [--region X] [--profile Y] s3 <cmd> [--only-show-errors] ...
RUN printf '%s\n' '#!/bin/bash' \
    '# Parse global AWS options before s3 subcommand' \
    'while [[ "$1" == --* ]] && [[ "$1" != "s3" ]]; do' \
    '    case "$1" in' \
    '        --region) export AWS_REGION="$2"; shift 2 ;;' \
    '        --profile) export AWS_PROFILE="$2"; shift 2 ;;' \
    '        *) shift ;;  # skip unknown global options' \
    '    esac' \
    'done' \
    'if [[ "$1" == "s3" ]]; then' \
    '    shift' \
    '    # Remove aws-cli-specific flags that s5cmd does not support' \
    '    args=()' \
    '    for arg in "$@"; do' \
    '        [[ "$arg" != "--only-show-errors" ]] && args+=("$arg")' \
    '    done' \
    '    # s5cmd uses "cat" for stdout output, not "cp source -"' \
    '    if [[ "${args[0]}" == "cp" ]] && [[ "${args[-1]}" == "-" ]]; then' \
    '        unset "args[-1]"  # remove trailing -' \
    '        args[0]="cat"     # change cp to cat' \
    '    fi' \
    '    exec /usr/local/bin/s5cmd "${args[@]}"' \
    'fi' \
    'echo "Unsupported aws command: $*" >&2' \
    'exit 1' > /usr/local/bin/aws && chmod +x /usr/local/bin/aws

# Copy only the cleaned conda environment
COPY --from=builder /opt/conda /opt/conda

WORKDIR /data
ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["build-fragments-h5", "--help"]
