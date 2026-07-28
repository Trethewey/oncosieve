# =============================================================================
# OncoSieve -- Docker Container
#
# Build:
#   docker build -t oncosieve:1.0 -f packaging/Dockerfile .
#
# Run:
#   docker run -v /path/to/data:/app/data -v /path/to/output:/app/output oncosieve:1.0
#
# Or use docker-compose:
#   docker-compose -f packaging/docker-compose.yml up
# =============================================================================

FROM ubuntu:22.04

LABEL maintainer="Dr Christopher Trethewey <christopher.trethewey@nhs.net>"
LABEL description="OncoSieve -- Pan-Cancer Somatic Variant Whitelist Curation"

# Avoid interactive prompts during build
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    bcftools \
    tabix \
    samtools \
    wget \
    gzip \
    && rm -rf /var/lib/apt/lists/*

# Verify versions
RUN python3 --version && bcftools --version | head -1 && tabix --version 2>&1 | head -1

# Set up application directory
WORKDIR /app

# Copy requirements first (Docker layer caching)
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy application code
COPY build_whitelist.py .
COPY post_process_whitelist.py .
COPY mutect2_rescue.py .
COPY pre_check.py .
COPY test_pipeline.py .
COPY config.yaml .
COPY settings.yaml.example .
COPY run_oncosieve.sh .
COPY parsers/ parsers/
COPY tools/ tools/
COPY panels/ panels/
COPY assets/ assets/

# Create output and intermediate directories
RUN mkdir -p output intermediate logs

# Copy settings.yaml if it exists, otherwise use example
RUN if [ ! -f settings.yaml ]; then cp settings.yaml.example settings.yaml; fi

# Make scripts executable
RUN chmod +x run_oncosieve.sh

# Data directory is expected to be mounted as a volume
VOLUME ["/app/data", "/app/output"]

# Default entrypoint runs the full pipeline
ENTRYPOINT ["bash", "run_oncosieve.sh"]

# Default: no arguments (use mounted data directory)
CMD []
