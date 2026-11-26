# Build environment of the prototype: Ubuntu 24.04 with GCC, GMP 6.3.0 and
# Python 3 from the Ubuntu archive, OpenSSL 3.5.9 (for ML-DSA-44) and RELIC
# 0.7.0 built from pinned commits.
#   docker build -t ect .
#   docker run --rm -it ect                     # a shell in src/
#   docker run --rm ect make tests              # the correctness tests
#   docker run --rm -v "$PWD:/ect" ect make benchmarks   # results written to the mounted repository
FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential cmake git ca-certificates perl libgmp-dev python3 util-linux \
 && rm -rf /var/lib/apt/lists/*

# OpenSSL 3.5.9 (tag openssl-3.5.9), installed in /opt/openssl.
ARG OPENSSL_COMMIT=45e844fa2a14ec92d146bd8f5778ac130b6625fb
RUN git clone https://github.com/openssl/openssl.git /tmp/openssl \
 && cd /tmp/openssl && git checkout -q $OPENSSL_COMMIT \
 && ./Configure --prefix=/opt/openssl --libdir=lib no-tests \
 && make -j"$(nproc)" && make install_sw \
 && rm -rf /tmp/openssl
ENV CPATH=/opt/openssl/include \
    LIBRARY_PATH=/opt/openssl/lib \
    LD_LIBRARY_PATH=/opt/openssl/lib \
    PATH=/opt/openssl/bin:$PATH

# RELIC 0.7.0: the last commit of the 0.7.0 series, with the options of the README.
ARG RELIC_COMMIT=9e6833250b476d79a4d985f176b6150553411ca1
RUN git clone https://github.com/relic-toolkit/relic.git /tmp/relic \
 && cd /tmp/relic && git checkout -q $RELIC_COMMIT \
 && mkdir build && cd build \
 && cmake -DALLOC=DYNAMIC -DBN_PRECI=8192 -DTESTS=0 -DBENCH=0 .. \
 && make -j"$(nproc)" && make install && ldconfig \
 && rm -rf /tmp/relic

COPY . /ect
WORKDIR /ect/src
RUN make clean && make && make test_dj test_dj_board test_rlwe test_eg_board bench_dj bench_rlwe bench_eg_board bench_trace micro_bench
CMD ["/bin/bash"]
