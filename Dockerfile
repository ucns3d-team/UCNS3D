# Start from the official Focal Fossa (20.04 LTS) image
ARG BASE_IMAGE=ubuntu:20.04
FROM ${BASE_IMAGE}

RUN apt-get update && \
    export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y \
    ssh \
    # wget \
    rpm \
    build-essential\
    zlib1g-dev \
    binutils-dev \
    gcc\
    g++\
    gfortran mpich \
    libblas-dev \
    liblapack-dev \
    autoconf \
    autotools-dev \
    # gnuplot \
    libfl-dev \
    libreadline-dev \
    libgmp-dev \
    libmpfr-dev \
    libfftw3-dev \
    libboost-system-dev \
    libboost-thread-dev \
    libcgal-dev \
    libmpc-dev \
    # ping \
    libstdc++5 \
    libiberty-dev \
    software-properties-common ;\
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y \
    git \
    cmake \
    mpich \
    libopenblas-dev ;\
    rm -rf /var/lib/apt/lists/*

WORKDIR /root/

#RUN mkdir -p /root/local ;\
#   cd /root/local \
#    cd /root/ ;\
#  pwd ;\
#    git clone https://github.com/KarypisLab/GKlib.git

#RUN make -C /root/GKlib prefix=/root/local config
#RUN	cd /root/GKlib && make install

#RUN cd /root/ ;\
#    git clone https://github.com/KarypisLab/METIS.git
#RUN make -C /root/METIS config shared=1 cc=gcc prefix=/root/local
#RUN	cd /root/METIS && make install


#RUN cd /root/ ;\
#    git clone https://github.com/KarypisLab/ParMETIS.git
#RUN make -C /root/ParMETIS config cc=mpicc prefix=/root/local
#RUN	cd /root/ParMETIS && make install
ADD src /root/CODE
RUN cd /root/CODE && make -f Makefile_docker clean all
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none

# Make & set a rundir & copy executable 
RUN mkdir -p /ucns3d_run
RUN cp /root/CODE/ucns3d_p /ucns3d_run

# Clean ubuntu dependencies
RUN apt-get clean && rm -rf /var/lib/apt/lists/*
# Remove src and metis dependencies
RUN rm -rf /root/CODE && rm -rf /root/GKlib && rm -rf /root/METIS && rm -rf /ParMETIS

WORKDIR /ucns3d_run
# Copy taylor green problem and run script
# TODO abstract testing out from dockerfile.
COPY tests/taylor_green_vortex/* /ucns3d_run/
COPY tests/execute-tests.sh /ucns3d_run/
RUN chmod +x execute-tests.sh
