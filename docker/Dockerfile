FROM conda/miniconda3@sha256:2141bfc4b60cc5eb332402207c80d884daa72fcefa48b50f3ccadd934d1f3d03

MAINTAINER loipf

#RUN apt-get update 
#RUN conda update --all

RUN apt-get update && apt install -y procps gzip curl

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge


### bioinformatic tools
RUN conda install \
		multiqc==1.9 \
		fastqc==0.11.9-0 \
		cutadapt==2.10 \
		kallisto==0.46.2 \
		r-base==3.6.3


RUN R -e 'install.packages("pacman",repos = "http://cran.us.r-project.org"); pacman::p_load(knitr, data.table, plotly, ggplot2, pheatmap, vsn, hexbin, LSD, cowplot, rmarkdown, tximport, DESeq2, rhdf5, dplyr, tidyr, rjson)'


