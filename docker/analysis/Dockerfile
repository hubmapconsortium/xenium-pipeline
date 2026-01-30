FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
 && apt-get -y install \
    autoconf \
    automake \
    bison \
    fastqc \
    flex \
    git \
    libfreetype6-dev \
    libxml2-dev \
    libtool \
    m4 \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    python-is-python3 \
    zlib1g-dev \
 && rm -rf /var/cache/apt/*

RUN python3 -m venv /opt/venv
ENV PATH=/opt/venv/bin:${PATH}

## Python packages
COPY docker/analysis/requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

COPY salmon-rnaseq/data /opt/data
COPY bin/common /opt
COPY bin/analysis /opt
COPY salmon-rnaseq/bin/analysis/plot_utils.py /opt

WORKDIR /opt
