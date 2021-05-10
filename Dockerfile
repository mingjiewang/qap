FROM continuumio/miniconda:4.3.14

LABEL maintainer="Mingjie_Wang huzai@sjtu.edu.cn"

ADD . /opt/qap

RUN conda config --add channels bioconda 
RUN conda config --set show_channel_urls yes 
RUN conda install cutadapt 
RUN apt-get install -y vim
RUN echo "alias ls='ls --color' " >> ~/.bashrc 
RUN echo "alias ll='ls -lah --color' " >> ~/.bashrc 
RUN echo "PS1='\[\\033[01;32m\]\u\[\\033[01;34m\]@\h:\w\$\[\\033[00m\]'" >> ~/.bashrc 
RUN echo 'export PATH=/opt/qap:$PATH' >> ~/.bashrc 
RUN cd /opt/qap/bin/3rdPartyTools/ 
RUN echo 'export JAVA_HOME=/opt/qap/bin/3rdPartyTools/jdk' >> ~/.bashrc 
RUN echo 'export JRE_HOME=/opt/qap/bin/3rdPartyTools/jdk/jre' >> ~/.bashrc 
RUN echo "export CLASSPATH=.:/opt/qap/bin/3rdPartyTools/jdk/lib/dt.jar:/opt/qap/bin/3rdPartyTools/jdk/lib/tools.jar" >> ~/.bashrc
RUN echo "cat /opt/qap/lib/bashrcFile >> ~/.bashrc"
RUN . ~/.bashrc  
RUN ln -sf /opt/qap/bin/3rdPartyTools/jdk/bin/* /usr/bin/ 
RUN pip install numpy && pip install biopython weblogo && apt install libgsl0ldbl gsl-bin 
RUN echo "deb http://cloud.r-project.org/bin/linux/debian jessie-cran34/" >> /etc/apt/sources.list 
RUN apt update 
RUN apt install -y --force-yes libcurl4-openssl-dev libssl-dev vim libatlas3-base r-base r-base-dev bioperl make gcc autoconf automake g++ pkg-config libpcre3 libpcre3-dev libgcrypt11-dev zlib1g-dev gsl-bin libgsl0-dev libgsl0ldbl 
RUN cd /opt/qap 
RUN R CMD javareconf 
RUN Rscript /opt/qap/bin/Rscripts/install_r_pkg.R 
RUN git clone https://github.com/vqv/ggbiplot 
RUN R CMD INSTALL ggbiplot && rm -rf ggbiplot 
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
RUN cp /opt/qap/README.md ~/
RUN cd ~/

