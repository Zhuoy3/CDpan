# Description:
# Author: Zhuo Yue
# Date: 2021-04-20
# Version 0.0.1

FROM centos:latest

LABEL maintainer="zhuoy"

USER root

WORKDIR /opt/

# Install dependencies of Python3
RUN yum -y install zlib zlib-devel
RUN yum -y install bzip2 bzip2-devel
RUN yum -y install ncurses ncurses-devel
RUN yum -y install readline readline-devel
RUN yum -y install openssl openssl-devel
RUN yum -y install xz xz-devel
RUN yum -y install sqlite sqlite-devel
RUN yum -y install gdbm gdbm-devel
RUN yum -y install tk tk-devel
RUN yum -y install gcc
RUN yum -y install make

# Install Python3 (needed by cutadapt)
COPY Python-3.9.4.tgz ./
RUN tar zxf Python-3.9.4.tgz
RUN rm -rf Python-3.9.4.tgz
WORKDIR /opt/Python-3.9.4/
RUN ./configure --with-ssl --prefix=/service/python3
RUN make
RUN make install
ENV PATH "$PATH:/service/python3/bin/"
WORKDIR /opt/
RUN rm -rf /opt/Python-3.9.4/

# Install cutadapt
RUN python3 -m pip install --user --upgrade cutadapt
ENV PATH "$PATH:/root/.local/bin"
RUN echo $PATH

# Install java (needed by FastQC)
RUN yum -y install java

# Install FastQC
COPY fastqc_v0.11.9.zip ./
RUN yum -y install zip
RUN unzip fastqc_v0.11.9.zip
RUN chmod 755 /opt/FastQC/fastqc
RUN ln -s /opt/FastQC/fastqc /usr/bin/fastqc
RUN rm -rf fastqc_v0.11.9.zip
RUN yum -y install perl

# Install TrimGalore
COPY TrimGalore-0.6.6.tar.gz ./
RUN tar xvzf TrimGalore-0.6.6.tar.gz
RUN ln -s /opt/TrimGalore-0.6.6/trim_galore /usr/bin/trim_galore
RUN rm -rf TrimGalore-0.6.6.tar.gz
RUN export LC_ALL=C
