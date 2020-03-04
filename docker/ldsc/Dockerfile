FROM continuumio/anaconda

WORKDIR /work

RUN apt update && apt -y install g++

RUN git clone https://github.com/bulik/ldsc.git &&  \
    cd ldsc && \
    conda env create --file environment.yml

