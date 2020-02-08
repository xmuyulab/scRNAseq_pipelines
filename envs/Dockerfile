FROM ubuntu:16.04
LABEL authors="wangshun@aginome.com" \
    description="Docker base image for scRNA pipeline by YuLab@XMU"

#ADD sources.list /etc/apt

RUN apt-get clean && \
    apt-get update --fix-missing && \
    apt-get install -y --fix-missing \
        build-essential \
        zlib1g-dev unzip \
        libncurses5-dev \
        curl \
        wget \
        python \
        python-pip \
        python-dev \
        cmake \
        libboost-all-dev \
        perl \
        perl-base \
        libxml2-dev \
        libcurl4-gnutls-dev \
        apt-transport-https \
        default-jre \
        default-jdk \
        less \
        vim \
        libtbb-dev \
        git \
        bzip2 \
        gcc \
        g++ \
        make \
        software-properties-common \
        tzdata \
        openssh-server \
        openssh-client

RUN cd /usr/local/bin && wget -qO- https://get.nextflow.io | bash && \
        chmod 755 ./nextflow && \
        nextflow -h && cp -r /root/.nextflow /.nextflow && \
        chmod -R 777 /.nextflow
RUN cd /opt &&\
        wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3.sh && \
        bash ./Miniconda3.sh -b -p /opt/miniconda3
ENV PATH $PATH:/opt/miniconda3/bin:/opt/Drop-seq_tools-2.3.0
ADD . /opt
RUN cd /opt && conda env create -f /opt/cellranger.yml && \
        tar zxf cellranger-3.0.2.tar.gz && rm -f cellranger-3.0.2.tar.gz && \
        mv -v /opt/cellranger-3.0.2/* /opt/miniconda3/envs/cellranger/bin/. && \
        rm -rf /opt/cellranger-3.0.2
RUN conda env create -f /opt/dropEst.yml
RUN cd /opt && conda env create -f /opt/Dropseq.yml && \
        unzip Drop-seq_tools-2.3.0.zip && rm -f Drop-seq_tools-2.3.0.zip
RUN conda env create -f /opt/scPipe.yml
RUN conda env create -f /opt/umis.yml
RUN conda env create -f /opt/UMI-tools.yml
RUN cd /opt && conda env create -f /opt/zUMIs.yml && \
        tar zxf zUMIs-2.5.4.tar.gz && rm -f zUMIs-2.5.4.tar.gz && \
        mv -v /opt/zUMIs-2.5.4/* /opt/miniconda3/envs/zUMIs/bin/. && \
        rm -rf /opt/zUMIs-2.5.4 && \
        ln -s /opt/miniconda3/envs/zUMIs/lib/R/modules/lapack.so /opt/miniconda3/envs/zUMIs/lib/libRlapack.so
RUN rm -rf  /opt/miniconda3/pkgs
