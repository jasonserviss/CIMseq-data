FROM rocker/tidyverse:3.5.2

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

RUN Rscript -e "source('https://raw.githubusercontent.com/jasonserviss/install/master/install_cran.R');install_cran(c('openxlsx/4.1.0', 'googledrive/0.1.3', 'lubridate/1.7.4', 'readxl/1.1.0'))"

# Clone and install EngeMetadata
RUN git clone https://github.com/EngeLab/EngeMetadata.git /home/EngeMetadata
RUN Rscript -e "devtools::install('/home/EngeMetadata', dependencies = FALSE)"

# Clone and install CIMseq.data
RUN touch /tmp.txt
RUN git clone https://github.com/jasonserviss/CIMseq.data.git /home/CIMseq.data
RUN Rscript -e "devtools::install('/home/CIMseq.data', dependencies = FALSE)"

# Run data scripts
WORKDIR /home/CIMseq.data

