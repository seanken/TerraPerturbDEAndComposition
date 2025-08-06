FROM rocker/tidyverse:4.4
WORKDIR /app

COPY scripts/GetDE.Gene.R /app/GetDE.Gene.R 
COPY scripts/GetDE.guide.R  /app/GetDE.guide.R
COPY scripts/TestClusterProportions.R /app/TestClusterProportions.R 

RUN sudo apt-get update && \
    apt-get install -y libglpk-dev libhdf5-dev hdf5-tools && \
    rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('remotes', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('aod',repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('qs',repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('SingleCellExperiment')"
RUN R -e "install.packages('nebula')"
RUN R -e "BiocManager::install('glmGamPoi')"
RUN R -e "remotes::install_github('bnprks/BPCells/r')"
RUN R -e "library(aod);library(BPCells);library(nebula);library(Seurat);library(glmGamPoi)"
