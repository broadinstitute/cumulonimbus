# Ubuntu with tabix, plink, python3, chowser and ecaviar

FROM ubuntu:18.04

WORKDIR /work

RUN apt update && apt -y install tabix

RUN apt update && apt -y install wget unzip

# install PLINK v1.90b6.7 64-bit (2 Dec 2018)
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20181202.zip && \
	unzip plink_linux_x86_64_20181202.zip && \
	chmod +x plink && \
    ln -s /work/plink /usr/local/bin/plink

RUN apt update && apt -y install python3

# install eCAVIAR prerequisites
RUN apt update && apt -y install g++ cmake git libgsl-dev liblapack3 liblapack-dev libblas3 libblas-dev

# install eCAVIAR
RUN git clone https://github.com/fhormoz/caviar.git
RUN cd caviar/CAVIAR-C++ && make && chmod +x eCAVIAR

# link to eCAVIAR
RUN ln -s /work/caviar/CAVIAR-C++/eCAVIAR /usr/local/bin/eCAVIAR

# install JVM
RUN apt update && apt -y install openjdk-11-jdk-headless

# install SBT and lintian
RUN apt -y update && \
    apt -y install gnupg2 fakeroot lintian && \
    echo "deb https://dl.bintray.com/sbt/debian /" | tee -a /etc/apt/sources.list.d/sbt.list && \
    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 2EE0EA64E40A89B84B2DF73499E82A75642AC823 && \
    apt-get update && \
    apt-get install sbt

# install Chowser 1.6.2
RUN git clone https://github.com/broadinstitute/chowser.git && \
    cd chowser && \
    sbt debian:lintian && \
    apt -y install ./target/chowser_1.6.2_all.deb
