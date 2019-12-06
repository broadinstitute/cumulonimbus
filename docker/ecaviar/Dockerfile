# Ubuntu with tabix, plink, python3, chowser and ecaviar

FROM ubuntu:18.04

WORKDIR /work

RUN apt update && apt -y install tabix

RUN apt update && apt -y install wget unzip

# install PLINK v1.90b6.7 64-bit (2 Dec 2018)
RUN mkdir tmp_plink && \
    cd tmp_plink && \
    wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20181202.zip && \
    unzip plink_linux_x86_64_20181202.zip && \
    cd .. && \
    mv tmp_plink/plink plink && \
    rm -r tmp_plink && \
	  chmod +x plink && \
    ln -s /work/plink /usr/local/bin/plink

RUN apt update && apt -y install python3

# install eCAVIAR
RUN apt update && apt -y install g++ cmake git libgsl-dev liblapack3 liblapack-dev libblas3 libblas-dev && \
    git clone https://github.com/fhormoz/caviar.git && \
    cd caviar/CAVIAR-C++ && make && chmod +x eCAVIAR && \
    cd /work && \
    mv /work/caviar/CAVIAR-C++/eCAVIAR eCAVIAR && \
    chmod u+x eCAVIAR && \
    rm -r caviar && \
    ln -s /work/eCAVIAR /usr/local/bin/eCAVIAR && \
    apt -y remove cmake g++ && \
    apt -y autoremove

# install JVM
RUN apt update && apt -y install openjdk-11-jdk-headless

# install SBT, lintian, chowser 1.7.1
RUN apt -y update && \
    apt -y install gnupg2 fakeroot lintian && \
    echo "deb https://dl.bintray.com/sbt/debian /" | tee -a /etc/apt/sources.list.d/sbt.list && \
    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 2EE0EA64E40A89B84B2DF73499E82A75642AC823 && \
    apt-get update && \
    apt-get install sbt && \
    git clone https://github.com/broadinstitute/chowser.git && \
    cd chowser && \
# workaround for SBT issue https://github.com/sbt/sbt/issues/5270
    sbt sbtVersion && \
    sbt debian:lintian && \
    apt -y install ./target/chowser_1.7.1_all.deb && \
    cd /work && \
    rm -r chowser && \
    apt -y remove sbt gnupg2 fakeroot lintian && \
    apt -y autoremove


    
