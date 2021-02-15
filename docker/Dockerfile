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
		deeptools==3.5.0-0 \
		cutadapt==2.10 \
		kallisto==0.46.2


## R install 
## perl 


 apt-get -y install r-base
		
#
# Install Sleuth 
# https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/
#
RUN R -e 'source("http://bioconductor.org/biocLite.R"); library(BiocInstaller); biocLite(c("XML","biomaRt")); biocLite("rhdf5"); install.packages("devtools", repos="http://cloud.r-project.org/"); devtools::install_github("pachterlab/sleuth")'




