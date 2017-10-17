FROM continuumio/miniconda:4.3.14

MAINTAINER Jianfeng Li lee_jianfeng@sjtu.edu.cn

ADD . /opt/qap

ENV JAVA_HOME /opt/qap/bin/3rdPartyTools/jdk1.8.0_131

RUN conda config --add channels bioconda && \
    conda config --set show_channel_urls yes && \
    conda install cutadapt && \
    cd /opt/qap/bin/3rdPartyTools/ && \
    tar -xzvf jdk1.8.0_131.tar.gz && \
    echo "export JAVA_HOME=/opt/qap/bin/3rdPartyTools/jdk1.8.0_131" >> ~/.bashrc && \
    echo "export JRE_HOME=$JAVA_HOME/jre" >> ~/.bashrc && \
    echo "export CLASSPATH=.:$JAVA_HOME/lib/dt.jar:$JAVA_HOME/lib/tools.jar" >> ~/.bashrc && \
    ln -s $JAVA_HOME/bin/* /usr/bin/ &&\
    . ~/.bashrc && \
    pip install numpy && pip install biopython weblogo && apt install libgsl0ldbl gsl-bin &&\
    echo "deb http://cloud.r-project.org/bin/linux/debian jessie-cran34/" >> /etc/apt/sources.list && \
    apt update && \
    apt install -y --force-yes libatlas3-base r-base r-base-dev bioperl make gcc autoconf automake g++ pkg-config \
    libpcre3 libpcre3-dev libgcrypt11-dev zlib1g-dev && \
    cd /opt/qap && \
    R CMD javareconf && \
    Rscript bin/3rdPartyTools/install_r_pkg.R &&\
    git clone https://github.com/vqv/ggbiplot && \
    R CMD INSTALL ggbiplot && rm -rf ggbiplot &&\
    cpan -fi Perl4::CoreLibs
