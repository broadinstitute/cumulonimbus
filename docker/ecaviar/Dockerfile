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

# install chowser 1.3.0
RUN wget https://personal.broadinstitute.org/oliverr/software/chowser/chowser_1.3.0_all.deb && \
    apt -y install ./chowser_1.3.0_all.deb
