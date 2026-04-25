FROM bioconductor/bioconductor_docker:devel
RUN R -e 'install.packages("remotes")'
RUN mkdir /pkg
COPY DESCRIPTION /pkg/DESCRIPTION
RUN R -e 'remotes::install_deps("/pkg", dependencies = TRUE)'
RUN R -e 'BiocManager::install("BiocCheck")'
