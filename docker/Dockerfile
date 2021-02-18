FROM continuumio/miniconda3@sha256:7838d0ce65783b0d944c19d193e2e6232196bada9e5f3762dc7a9f07dc271179

MAINTAINER loipf

RUN conda update --all

RUN apt-get update && apt install -y procps gzip curl chrpath

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge


### bioinformatic tools
RUN conda install \
		multiqc==1.9 \
		fastqc==0.11.9-0 \
		cutadapt==2.10 \
		kallisto==0.46.2

RUN apt install -y r-base
RUN apt install -y build-essential aptitude libxml2-dev libssl-dev libcurl4-openssl-dev libxtst6 libgl1-mesa-glx

RUN R -e 'install.packages("pacman",repos = "http://cran.us.r-project.org")'
RUN R -e 'pacman::p_load(BiocManager, knitr, data.table, plotly, ggplot2, pheatmap, hexbin, LSD, cowplot, rmarkdown, dplyr, tidyr, tximport,rhdf5, vsn, rjson)'

