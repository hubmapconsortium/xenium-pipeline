FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

ARG SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True

RUN apt-get update \
 && apt-get -y install \
    autoconf \
    automake \
    bison \
    fastqc \
    flex \
    git \
    libfreetype6-dev \
#    libigraph0-dev \
    libxml2-dev \
    libtool \
    m4 \
    pkg-config \
    python3-dev \
    python3-pip \
    python3-venv \
    python-is-python3 \
    zlib1g-dev \
 && rm -rf /var/cache/apt/*


RUN python3 -m venv /opt/venv
ENV PATH=/opt/venv/bin:${PATH}

## Python packages
COPY docker/conversion/requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

COPY bin/common /opt
COPY bin/conversion /opt

WORKDIR /opt

