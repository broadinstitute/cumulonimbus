# Ubuntu with tabix, plink, python3 and ecaviar

FROM ubuntu:18.04

WORKDIR /work

RUN apt update && apt -y install tabix

RUN apt update && apt -y install wget unzip

RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20181202.zip && \
	unzip plink_linux_x86_64_20181202.zip && \
	chmod +x plink && \
    ln -s /work/plink /usr/local/bin/plink

RUN apt update && apt -y install python3

# install cmake
RUN apt update && apt -y install g++ cmake git libgsl-dev liblapack3 libblas3

# install eCAVIAR
RUN git clone https://github.com/fhormoz/caviar.git
RUN cd caviar/CAVIAR-C++ && make && chmod +x eCAVIAR

# link to eCAVIAR
RUN ln -s caviar/CAVIAR-C++/eCAVIAR /usr/local/bin/eCAVIAR
