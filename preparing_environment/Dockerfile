# Base Image
FROM biocontainers/biocontainers:latest
USER root
RUN apt-get update
RUN apt-get install docker.io -y
RUN apt-get install python-psutil -y
RUN apt-get install python-numpy -y
RUN apt-get install python-pandas -y
RUN apt-get install tabix
RUN conda install bcftools
RUN git clone https://github.com/mrinalmanu/gatk4_exome_scripts /home/NGS_pipeline/
WORKDIR "/home/NGS_pipeline"
