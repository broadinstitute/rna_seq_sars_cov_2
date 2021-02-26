FROM continuumio/miniconda3:4.7.12

SHELL ["/bin/bash", "-c"]
RUN mkdir /usr/share/man/man1/
# Download cellranger-5.0.0.tar.gz 1st
RUN apt-get update && apt-get install -y  \
    autoconf \
    build-essential \
    bzip2 \
    default-jre \
    git \
    libbz2-dev \
    libdb-dev \
    liblzma-dev \
    perl \
    unzip \
    wget \
    zlib1g \
    zlib1g-dev \
    zlibc

RUN pip install --upgrade pip && \
    pip install pandas igv-reports

ADD https://github.com/deweylab/RSEM/archive/v1.3.3.zip /software/
RUN unzip /software/v1.3.3.zip -d /software/ && cd /software/RSEM-1.3.3/ && make && make install

RUN wget https://github.com/alexdobin/STAR/archive/2.7.2b.tar.gz && \
    tar xvf 2.7.2b.tar.gz -C /software && \
    rm 2.7.2b.tar.gz

ENV SAMTOOLS_VERSION=1.9
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 -C /software && \
    cd /software/samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
    cd ../ && ./configure --prefix=/software --without-curses && make && make install

ENV GATK_VERSION=4.1.9.0
RUN wget https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    unzip gatk-${GATK_VERSION}.zip -d /software/ && \
    rm gatk-${GATK_VERSION}.zip

RUN git clone --recursive https://github.com/trinityrnaseq/bamsifter.git && \
    mv bamsifter /software/bamsifter && cd /software/bamsifter && make

RUN wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 && \
    tar xf minimap2-2.17_x64-linux.tar.bz2 && \
    rm minimap2-2.17_x64-linux.tar.bz2 && \
    mv minimap2-2.17_x64-linux /software/minimap2-2.17_x64-linux
ADD cellranger-5.0.0.tar.gz /software
RUN apt-get -qq -y remove autoconf build-essential dpkg-dev git lsb-release unzip && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

ENV PATH=/software/minimap2-2.17_x64-linux/:/software/gatk-${GATK_VERSION}:/software/cellranger-5.0.0:/software/samtools-${SAMTOOLS_VERSION}/:/software/bin:/software/RSEM-1.3.3:/software/STAR-2.7.2b/bin/Linux_x86_64_static/:/software/bamsifter:$PATH:

