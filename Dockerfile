FROM rocker/tidyverse:3.4.3

MAINTAINER Jason Serviss <jason.serviss@ki.se>

# System dependencies for required R packages
RUN  rm -f /var/lib/dpkg/available \
  && rm -rf  /var/cache/apt/* \
  && apt-get update -qq \
  && apt-get install -y --no-install-recommends \
    ca-certificates \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    git

# R dependencies
RUN Rscript -e "install.packages(c('devtools','knitr','rmarkdown','shiny','RCurl'), repos = 'https://cran.rstudio.com')"

RUN Rscript -e "source('https://cdn.rawgit.com/road2stat/liftrlib/aa132a2d/install_cran.R');install_cran(c('openxlsx/4.0.17', 'googledrive/0.1.1'))"

# Clone and install sp.scRNAseqData
RUN git clone https://github.com/jasonserviss/sp.scRNAseqData.git /home/sp.scRNAseqData

RUN Rscript -e "devtools::install('/home/sp.scRNAseqData')"

# Run data scripts
RUN Rscript -e "source('/home/sp.scRNAseqData/inst/rawData/processRaw.R')"

WORKDIR /home
