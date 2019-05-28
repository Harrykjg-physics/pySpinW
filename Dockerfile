# Download and install Matlab Compiler Runtime v9.6 (2019a)
#
# This docker file will configure an environment into which the Matlab compiler
# runtime will be installed and in which stand-alone matlab routines (such as
# those created with Matlab's deploytool) can be executed.
#
# See http://www.mathworks.com/products/compiler/mcr/ for more info.

FROM debian:stretch-slim

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -q update \
    && apt-get install -q -y --no-install-recommends \
         xorg \
         unzip \
         wget \
         curl \
         bzip2 \
         ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install conda
RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local \
    && rm -rf /tmp/miniconda.sh \
    && conda install -y python=3 \
    && conda update conda \
    && apt-get -qq -y remove curl bzip2 \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && conda clean --all --yes

ENV PATH /opt/conda/bin:$PATH
    
# Install the MCR dependencies and some things we'll need and download the MCR
# from Mathworks -silently install it
RUN mkdir /mcr-install \
    && mkdir /opt/mcr \
    && cd /mcr-install \
    && wget -q http://ssd.mathworks.com/supportfiles/downloads/R2019a/Release/2/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2019a_Update_2_glnxa64.zip \
    && unzip -q MATLAB_Runtime_R2019a_Update_2_glnxa64.zip \
    && rm -f MATLAB_Runtime_R2019a_Update_2_glnxa64.zip \
    && ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent \
    && cd / \
    && rm -rf mcr-install

# Configure environment variables for MCR
ENV LD_LIBRARY_PATH /opt/mcr/v96/runtime/glnxa64:/opt/mcr/v96/bin/glnxa64:/opt/mcr/v96/sys/os/glnxa64
ENV XAPPLRESDIR /opt/mcr/v96/X11/app-defaults

RUN conda create -y -n pySpinW python=3.7 anaconda \
    && conda init bash

# docker build -t spinw . && docker run -p 8888:8888 -v /Users/simonward/PycharmProjects/pySpinW2:/opt/notebooks -t -i spinw bash -c "source activate pySpinW && jupyter notebook --notebook-dir=/opt/notebooks --ip='0.0.0.0' --port=8888 --no-browser --allow-root"
# bash -c "source activate pySpinW && jupyter notebook --notebook-dir=/opt/notebooks --ip='0.0.0.0' --port=8888 --no-browser --allow-root"