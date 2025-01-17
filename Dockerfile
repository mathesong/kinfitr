# Use the Rocker RStudio image as base
FROM rocker/rstudio:latest

# Disable authentication
ENV AUTH=none

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /home/rstudio/kinfitr

# Set environment variable to indicate Docker environment
ENV DOCKER=true

# Copy project files
COPY . /home/rstudio/kinfitr/

# Make sure RStudio user owns the project files
RUN chown -R rstudio:rstudio /home/rstudio/kinfitr

# Install build dependencies first with parallel installation
RUN Rscript -e '\
    repos <- c(CRAN = "https://cloud.r-project.org"); \
    options(Ncpus = parallel::detectCores()); \
    install.packages(c("xml2", "roxygen2", "remotes", "devtools"), repos = repos, dependencies = TRUE)'

# Install R package dependencies from DESCRIPTION file with parallel installation
RUN Rscript -e '\
    repos <- c(CRAN = "https://cloud.r-project.org"); \
    options(Ncpus = parallel::detectCores()); \
    if (!require("remotes")) install.packages("remotes", repos = repos); \
    if (file.exists("DESCRIPTION")) { \
        deps <- read.dcf("DESCRIPTION", fields = c("Imports", "Depends")); \
        packages <- unlist(strsplit(paste(deps[!is.na(deps)], collapse = ","), ",")); \
        packages <- gsub("\\\\s+", "", packages); \
        packages <- packages[packages != "R"]; \
        install.packages(packages, repos = repos, dependencies = TRUE); \
        remotes::install_deps(dependencies = TRUE, upgrade = "never", Ncpus = parallel::detectCores()); \
    }'

# Expose port 8787 for RStudio
EXPOSE 8787