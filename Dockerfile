# syntax=docker/dockerfile:1.7
#
# mirna-kras — full pipeline container
#
# Builds everything from source at build time:
#   - MC-Flashfold (C → Linux ELF via gcc)   ← ./docker/build-context/MC-Flashfold-v37.0/
#   - RINexus (Java → .class via Maven)      ← ./docker/build-context/rimap_java_code/
#
# Quick start:
#
#   docker build -t mirna-kras .
#   docker run --rm -p 8000:8000 mirna-kras                              # mkdocs serve
#   docker run --rm -v $(pwd)/work:/app/work mirna-kras bash            # shell for the pipeline
#
# Both build-context bundles are committed to the repo, so the image is
# reproducible without any external download.

ARG PYTHON_VERSION=3.14-slim-bookworm
ARG JAVA_VERSION=21


# ──────────────────────────────────────────────────────────────────────────────
# Stage 1: compile MC-Flashfold from C source
# Source: ./docker/build-context/MC-Flashfold-v37.0/  (src/, tables/, tablesPermissives/)
# ──────────────────────────────────────────────────────────────────────────────
FROM debian:bookworm-slim AS flashfold-build
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        gcc \
        libc6-dev \
        make \
        ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /build/mcflashfold
COPY docker/build-context/MC-Flashfold-v37.0/ ./
RUN gcc -O3 src/flashfold.c -o mcff -lm \
 && gcc -O3 src/db2cm.c  -o db2cm -lm \
 && echo "MC-Flashfold compiled from C source" \
 && ls -lh mcff db2cm


# ──────────────────────────────────────────────────────────────────────────────
# Stage 2: compile RINexus Java classes with Maven
# Source: ./docker/build-context/rimap_java_code/  (pom.xml + src/)
#
# Maven dependencies (declared in the bundled pom.xml):
#   - spring-boot-starter-web
#   - spring-boot-starter-data-jpa
#   - jackson-databind
#   - commons-csv / commons-lang3 / commons-io / commons-codec
#   - lombok (provided)
#   - h2
#   - slf4j-api
# ──────────────────────────────────────────────────────────────────────────────
FROM debian:bookworm-slim AS rinexus-build
ARG JAVA_VERSION=21
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        openjdk-${JAVA_VERSION}-jdk-headless \
        maven \
        ca-certificates \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /build/rinexus
COPY docker/build-context/rimap_java_code/ ./
RUN mvn -B -q package -DskipTests \
 && cp -r target/classes /build/rinexus-classes \
 && find /build/rinexus-classes -name '*.class' | wc -l | xargs -I{} echo "RINexus: {} classes compiled"


# ──────────────────────────────────────────────────────────────────────────────
# Stage 3: runtime image — Python 3.14 + Java 21 + repo
# ──────────────────────────────────────────────────────────────────────────────
FROM python:${PYTHON_VERSION}
ARG JAVA_VERSION=21
ENV DEBIAN_FRONTEND=noninteractive

# Java 21 for RINexus, build tools for legacy scripts, libgomp for HDBSCAN's
# OpenMP runtime.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        openjdk-${JAVA_VERSION}-jdk-headless \
        git \
        make \
        gcc \
        g++ \
        libstdc++6 \
        libgomp1 \
        ca-certificates \
        wget \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Python deps — copy first so this layer caches.
COPY requirements.txt ./
RUN pip install --no-cache-dir --upgrade pip \
 && pip install --no-cache-dir -r requirements.txt


# MC-Flashfold binaries + energy tables (from Stage 1)
COPY --from=flashfold-build /build/mcflashfold/mcff               /usr/local/bin/mcff
COPY --from=flashfold-build /build/mcflashfold/db2cm             /usr/local/bin/db2cm
COPY --from=flashfold-build /build/mcflashfold/tables              /opt/MC-Flashfold/tables/
COPY --from=flashfold-build /build/mcflashfold/tablesPermissives  /opt/MC-Flashfold/tablesPermissives/
ENV MCFLASHFOLD_PATH=/opt/MC-Flashfold


# RINexus Java classes (from Stage 2)
COPY --from=rinexus-build /build/rinexus-classes/ /app/archive/rimap_java_code/target/classes/
ENV JAVA_CLASSPATH=/app/archive/rimap_java_code/target/classes


# Repo contents
COPY scripts/ ./scripts/
COPY docs/    ./docs/
COPY data/    ./data/
COPY mkdocs.yml README.md LICENSE ./


# Volume target for outputs (figures, intermediate databases)
VOLUME ["/app/work"]

# mkdocs dev server port
EXPOSE 8000

# Default command — overridable via `docker run mirna-kras <cmd>`
CMD ["mkdocs", "serve", "--dev-addr=0.0.0.0:8000"]
