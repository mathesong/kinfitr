Bootstrap: docker
From: rocker/rstudio:4.4.2

%post
    # System dependencies
    apt-get update && apt-get install -y \
        libxml2-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libgit2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        && rm -rf /var/lib/apt/lists/*

    # Install R packages
    R --quiet -e "install.packages(c('devtools', 'roxygen2', 'testthat', 'knitr', 'rmarkdown'), repos='https://cloud.r-project.org/')"
    
    # Clean up
    rm -rf /tmp/downloaded_packages

%environment
    export PATH=/usr/lib/rstudio-server/bin:$PATH

%startscript
    exec /init

%runscript
    exec /init 