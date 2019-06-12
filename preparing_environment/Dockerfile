# Base Image
FROM biocontainers/biocontainers:latest
USER root
RUN conda install bcftools=1.3.1
RUN apt-get update
RUN apt-get install docker.io -y
RUN apt-get install python-psutil -y
RUN apt-get install python-numpy -y
RUN apt-get install python-pandas -y
RUN apt-get install tabix
RUN git clone https://github.com/ptriska/NGS_pipeline/ /home/NGS_pipeline/
WORKDIR "/home/NGS_pipeline"
