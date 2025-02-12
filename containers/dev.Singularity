Bootstrap: docker
From: rocker/r-ver:4.4.2

%post
    # System updates and basic tools
    apt-get update && apt-get install -y \
        build-essential \
        git \
        curl \
        wget \
        r-base \
        r-base-dev \
        && rm -rf /var/lib/apt/lists/*

    # Install R packages
    R --quiet -e "install.packages(c('devtools', 'roxygen2', 'testthat'), repos='https://cloud.r-project.org/')"

%environment
    export LC_ALL=C

%runscript
    exec /bin/bash "$@" 