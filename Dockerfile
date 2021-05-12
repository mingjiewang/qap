FROM continuumio/miniconda:4.3.14

# Metadata
LABEL maintainer="Mingjie_Wang <huzai@sjtu.edu.cn>, XU_Nuo <xunuo_promise@outlook.com>"

# Copy project files
COPY . /opt/qap

# Setup conda 
RUN conda config --add channels bioconda 
RUN conda config --set show_channel_urls yes 
RUN conda install cutadapt

# Install tools required for project
RUN apt-get update && apt install -y --force-yes \
    dirmngr \
    apt-transport-https \
    ca-certificates \
    software-properties-common \
    gnupg2
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian jessie-cran34/'
RUN apt update && apt install -y --force-yes \
    autoconf \
    automake \
    bioperl \
    gcc \
    gsl-bin \
    g++ \
    libatlas3-base \
    libcurl4-openssl-dev \
    libgcrypt11-dev \
    libgsl0ldbl \
    libgsl0-dev \
    libpcre3 \
    libpcre3-dev \
    libssl-dev \
    make \
    pkg-config \
    r-base \
    r-base-dev \
    zlib1g-dev \
    vim \
    && rm -rf /var/lib/apt/lists/*

# Configure Bash
RUN echo "alias ls='ls --color' " >> /root/.bashrc 
RUN echo "alias ll='ls -lah --color' " >> /root/.bashrc 
RUN echo "PS1='\[\\033[01;32m\]\u\[\\033[01;34m\]@\h:\w\$\[\\033[00m\]'" >> /root/.bashrc 
RUN echo 'export PATH=/opt/qap:$PATH' >> /root/.bashrc 
WORKDIR /opt/qap/bin/3rdPartyTools/ 
RUN echo 'export JAVA_HOME=/opt/qap/bin/3rdPartyTools/jdk' >> /root/.bashrc 
RUN echo 'export JRE_HOME=/opt/qap/bin/3rdPartyTools/jdk/jre' >> /root/.bashrc 
RUN echo "export CLASSPATH=.:/opt/qap/bin/3rdPartyTools/jdk/lib/dt.jar:/opt/qap/bin/3rdPartyTools/jdk/lib/tools.jar" >> /root/.bashrc
RUN echo "cat /opt/qap/lib/bashrcFile >> /root/.bashrc"
RUN . /root/.bashrc

# Create soft link for Java
RUN chmod 755 /opt/qap/bin/3rdPartyTools/jdk/bin/*
RUN ln -sf /opt/qap/bin/3rdPartyTools/jdk/bin/* /usr/bin/

# Install Python dependency
RUN pip install numpy \
    && pip install biopython==1.76 \
    weblogo 

# Install R dependency
WORKDIR /opt/qap 
RUN R CMD javareconf 
RUN Rscript /opt/qap/bin/Rscripts/install_r_pkg.R 
RUN git clone https://github.com/vqv/ggbiplot 
RUN R CMD INSTALL ggbiplot \
    && rm -rf ggbiplot 

# Install Perl dependency
RUN cpan -fi Perl4::CoreLibs 
RUN cpan -fi AppConfig 
RUN cpan -fi Config::General 
RUN cpan -fi Font::TTF::Font 
RUN cpan -fi Math::Bezier 
RUN cpan -fi Math::VecStat 
RUN cpan -fi Math::Round 
RUN cpan -fi Params::Validate 
RUN cpan -fi Readonly 
RUN cpan -fi Set::IntSpan 
RUN cpan -fi Statistics::Basic 
RUN cpan -fi Text::Format

# Wrap up
RUN cp /opt/qap/README.md /root
WORKDIR /root

