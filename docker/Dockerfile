# syntax=docker/dockerfile:1.7
#
# AMFlow.cpp Docker image — bundles FLINT 3.4 + Kira 3.1 + FireFly +
# amflow_cli for out-of-the-box multi-loop Feynman integral computation.
#
# Fermat (closed-source, free academic use, NOT redistributable) is the
# only piece the user must provide at runtime via bind-mount; see
# docs/DOCKER.md.
#
# Build:
#   docker build -f docker/Dockerfile -t amflow-cpp:dev .
# Run:
#   docker run -v /path/to/Ferl7:/usr/share/Ferl7:ro -v "$(pwd)":/work \
#              amflow-cpp:dev \
#              /opt/amflow-cpp/tools/bench/cutbubble_1L_eps001_black_box_amflow_cpp.json \
#              /work/out.json

# ============================================================
# Stage 1 — build
# ============================================================
FROM ubuntu:22.04 AS build

ARG DEBIAN_FRONTEND=noninteractive
ARG FLINT_VERSION=3.4.0
ARG KIRA_VERSION=3.1

# Build-time deps:
#   build-essential / cmake / ninja / pkg-config / meson — toolchain
#   git / curl / ca-certificates                          — fetching sources
#   libgmp-dev / libmpfr-dev                              — arithmetic backends
#   libcln-dev / libginac-dev                             — Kira symbolic algebra
#   libyaml-cpp-dev                                       — Kira config parser
#   libjemalloc-dev / zlib1g-dev                          — Kira optional+core
#   libgtest-dev / nlohmann-json3-dev                     — amflow_cli build deps
#   python3                                               — Meson dependency
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake ninja-build pkg-config meson \
        autoconf automake libtool \
        git curl ca-certificates \
        libgmp-dev libmpfr-dev \
        libcln-dev libginac-dev \
        libyaml-cpp-dev libjemalloc-dev zlib1g-dev \
        libgtest-dev nlohmann-json3-dev \
        python3 \
    && rm -rf /var/lib/apt/lists/*
# Note: libkyotocabinet-dev intentionally omitted.  Kira 3.1's
# src/pyred/keyvaluedb.h includes `<kcpolydb.h>` inside an
# `#ifdef PYRED_KCDB` block placed BEFORE its `#include "pyred/config.h"`,
# so the macro isn't yet defined when the preprocessor decides whether
# to pull in the kyotocabinet header.  Kyoto is an optional pyRed
# backend; AMFlow.cpp doesn't depend on it, so we build with
# `-Dkyotocabinet=false` and save ~10MB of apt deps too.

# ---- FLINT 3.4 from source ----
# Ubuntu 22.04 ships FLINT 2.8.4; amflow_cli needs >= 3.4 for
# fmpz_mpoly_q_set_str_pretty and fmpz_mpoly_q_evaluate_acb.  Install to
# /opt/flint so it's easy to copy into the runtime stage cleanly.
RUN cd /tmp \
    && curl -fsSL --retry 8 --retry-delay 5 --retry-all-errors \
        "https://github.com/flintlib/flint/releases/download/v${FLINT_VERSION}/flint-${FLINT_VERSION}.tar.gz" \
        -o flint.tar.gz \
    && tar -xzf flint.tar.gz \
    && cd flint-${FLINT_VERSION} \
    && ./bootstrap.sh \
    && ./configure --prefix=/opt/flint --disable-static \
    && make -j"$(nproc)" \
    && make install \
    && cd / && rm -rf /tmp/flint-${FLINT_VERSION} /tmp/flint.tar.gz \
    # Expose FLINT to the dynamic linker (so kira / amflow_cli runtime
    # finds libflint.so.22) and to the build-time linker via
    # /etc/ld.so.conf.d (so meson's `find_library('flint')` from the
    # FireFly subproject succeeds; otherwise FireFly falls back to its
    # bundled flint2.wrap which tries to download FLINT 2.8.4 +
    # mesonpatches from kira.hepforge.org and breaks the build).
    && echo "/opt/flint/lib" > /etc/ld.so.conf.d/flint.conf \
    && ldconfig

ENV PKG_CONFIG_PATH=/opt/flint/lib/pkgconfig \
    LD_LIBRARY_PATH=/opt/flint/lib \
    LIBRARY_PATH=/opt/flint/lib \
    C_INCLUDE_PATH=/opt/flint/include \
    CPLUS_INCLUDE_PATH=/opt/flint/include \
    CMAKE_PREFIX_PATH=/opt/flint

# ---- Kira 3.1 + FireFly (subproject) ----
# Kira's meson auto-fetches FireFly from gitlab.com/firefly-library/firefly
# (revision kira-2) via subprojects/firefly.wrap.  We disable Kira's own
# flint2 fallback because we want to share the FLINT 3.4 install with
# amflow_cli — `flint=true` makes FireFly link the system FLINT (found via
# PKG_CONFIG_PATH=/opt/flint/lib/pkgconfig set above).
#
# Pre-generate `src/pyred/config.h` serially before the parallel ninja
# pass: Meson lists config.h as a source of the kira executable, but
# main.cpp's `#include "pyred/config.h"` dependency only takes effect
# via the .d file produced by the FIRST successful compile.  Under
# `-j$(nproc)`, ninja can therefore race — scheduling main.cpp.o before
# config.h's vcs_tag custom-target has run, so `VERSION` is undefined at
# preprocess time.  Forcing the config.h target first removes the race.
# Clone Kira AND pre-populate its FireFly subproject manually instead
# of letting Meson's wrap-git fetch it during `meson setup`.  Meson's
# wrap fetch has no retry logic, so a single gitlab.com TLS hiccup
# fails the whole build.  Cloning manually inside a retry loop is
# robust against transient network errors.
RUN cd /tmp \
    && for i in 1 2 3 4 5; do \
        git clone --branch ${KIRA_VERSION} https://gitlab.com/kira-pyred/kira.git && break \
        || { echo "kira clone attempt $i failed, retrying..."; sleep 5; rm -rf kira; }; \
       done \
    && [ -d kira ] || { echo "kira clone failed after 5 attempts"; exit 1; } \
    && cd kira/subprojects \
    && for i in 1 2 3 4 5; do \
        git clone --branch kira-2 https://gitlab.com/firefly-library/firefly.git && break \
        || { echo "firefly clone attempt $i failed, retrying..."; sleep 5; rm -rf firefly; }; \
       done \
    && [ -d firefly ] || { echo "firefly clone failed after 5 attempts"; exit 1; } \
    && cd /tmp/kira \
    && meson setup build \
        --prefix=/opt/kira \
        --buildtype=release \
        -Djemalloc=true \
        -Dkyotocabinet=false \
        -Dfirefly=true \
        -Dflint=true \
    && ninja -C build src/pyred/config.h \
    && ninja -C build -j"$(nproc)" \
    && ninja -C build install \
    && cd / && rm -rf /tmp/kira

# ---- amflow_cli ----
WORKDIR /src
COPY . /src
RUN cmake -S . -B build -G Ninja -DAMFLOW_BUILD_DRIVER=ON \
    && cmake --build build -j"$(nproc)"

# ============================================================
# Stage 2 — runtime
# ============================================================
FROM ubuntu:22.04 AS runtime

ARG DEBIAN_FRONTEND=noninteractive

# Runtime libraries only (no -dev packages).  Versions are pinned to
# Ubuntu 22.04's defaults so they match the build stage's link-time
# resolution; do not bump the base image without re-checking.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgmp10 libmpfr6 \
        libcln6 libginac11 \
        libyaml-cpp0.7 libjemalloc2 \
        zlib1g libstdc++6 \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# FLINT 3.4 runtime
COPY --from=build /opt/flint /opt/flint

# Kira install (bin/, lib/libfirefly*, share/, ...)
COPY --from=build /opt/kira /opt/kira

# amflow_cli binary lands at the path expected by run_json.cpp's
# kira_executable default — /usr/local/bin/kira is a symlink that wires
# up the "just works" experience where users don't need to override
# kira_executable in their bench JSON.
COPY --from=build /src/build/src/cli/amflow_cli /usr/local/bin/amflow_cli
RUN ln -s /opt/kira/bin/kira /usr/local/bin/kira

# Copy oracle benches + compare script (read-only inside image)
COPY --from=build /src/tools/bench /opt/amflow-cpp/tools/bench

# Library search paths.  Kira's meson.build forces libdir=lib (i.e. no
# debian multiarch x86_64-linux-gnu subdir) so /opt/kira/lib is the only
# Kira lib path we need.
ENV PATH=/opt/kira/bin:/usr/local/bin:/usr/bin:/bin \
    LD_LIBRARY_PATH=/opt/flint/lib:/opt/kira/lib

# Fermat must be mounted by the user (closed-source, non-redistributable).
# The hard-coded default in src/api/run_json.cpp is /usr/share/Ferl7/fer64;
# the entrypoint detects it and prints a clear error if missing.

# Entrypoint
COPY docker/entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

WORKDIR /work

ENTRYPOINT ["/entrypoint.sh"]
CMD []

LABEL org.opencontainers.image.title="AMFlow.cpp" \
      org.opencontainers.image.description="C++17 reimplementation of the AMFlow auxiliary-mass-flow algorithm for multi-loop Feynman integrals" \
      org.opencontainers.image.url="https://github.com/chang18/amflow-cpp" \
      org.opencontainers.image.source="https://github.com/chang18/amflow-cpp" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.documentation="https://github.com/chang18/amflow-cpp/blob/main/docs/DOCKER.md"
