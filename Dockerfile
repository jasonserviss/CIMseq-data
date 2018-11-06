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

RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R');install_cran(c('rlang/0.2.1', 'openxlsx/4.0.17', 'googledrive/0.1.1', 'dplyr/0.7.7', 'stringr/1.3.1', 'tibble/1.4.2', 'readr/1.1.1', 'purrr/0.2.5', 'lubridate/1.7.4', 'readxl/1.1.0', 'tidyr/0.8.1'))"

# Clone and install EngeMetadata
RUN git clone https://github.com/EngeLab/EngeMetadata.git /home/EngeMetadata
RUN Rscript -e "devtools::install('/home/EngeMetadata')"

# Clone and install sp.scRNAseqData
RUN git clone https://github.com/jasonserviss/sp.scRNAseqData.git /home/sp.scRNAseqData
RUN Rscript -e "devtools::install('/home/sp.scRNAseqData')"

# Run data scripts
WORKDIR /home/sp.scRNAseqData

