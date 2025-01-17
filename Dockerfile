# Use the Rocker RStudio image as base
FROM rocker/rstudio:latest

# Disable authentication and set default user credentials
ENV AUTH=none
ENV USER=rstudio
ENV PASSWORD=rstudio

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R package dependencies first
RUN Rscript -e '\
    repos <- c(CRAN = "https://cloud.r-project.org"); \
    options(Ncpus = parallel::detectCores()); \
    install.packages(c("xml2", "roxygen2", "remotes", "devtools", "usethis", "testthat", "pkgload"), \
        repos = repos, \
        dependencies = TRUE)'

# Set working directory
WORKDIR /home/rstudio/kinfitr

# Set environment variable to indicate Docker environment
ENV DOCKER=true

# Copy project files
COPY . /home/rstudio/kinfitr/

# Make sure RStudio user owns the project files
RUN chown -R rstudio:rstudio /home/rstudio/kinfitr

# Install R package dependencies from DESCRIPTION file
RUN Rscript -e '\
    repos <- c(CRAN = "https://cloud.r-project.org"); \
    options(Ncpus = parallel::detectCores()); \
    if (file.exists("DESCRIPTION")) { \
        deps <- read.dcf("DESCRIPTION", fields = c("Imports", "Depends")); \
        packages <- unlist(strsplit(paste(deps[!is.na(deps)], collapse = ","), ",")); \
        packages <- gsub("\\\\s+", "", packages); \
        packages <- packages[packages != "R"]; \
        install.packages(packages, repos = repos, dependencies = TRUE); \
    }'

# Expose port 8787 for RStudio
EXPOSE 8787