FROM rocker/r-ver:4.4.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages needed for development
RUN R -e 'install.packages(c("devtools", "roxygen2", "testthat", "lintr"))'

# Set working directory
WORKDIR /workspace 