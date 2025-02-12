# Use Rocker RStudio image as base
FROM rocker/rstudio:4.4.2

# Install RStudio Server
RUN /rocker_scripts/install_rstudio.sh

# Disable authentication and set default user credentials
ENV AUTH=none
ENV USER=rstudio
ENV PASSWORD=rstudio

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    # Additional system dependencies for R packages
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /home/rstudio/kinfitr

# Set environment variable to indicate Docker environment
ENV DOCKER=true

# Copy project files
COPY . /home/rstudio/kinfitr/

# Make sure RStudio user owns the project files
RUN chown -R rstudio:rstudio /home/rstudio/kinfitr

# Install devtools and then use it to install the package with all dependencies in parallel
RUN Rscript -e '\
    install.packages(c("devtools", "systemfonts", "textshaping", "ragg", "shiny", "miniUI", "pkgdown")); \
        #repos="https://cloud.r-project.org", \
        #dependencies=TRUE); \
    options(Ncpus = parallel::detectCores()); \
    devtools::install("/home/rstudio/kinfitr", dependencies=TRUE, upgrade="never")'

RUN Rscript -e 'cat("\nR Library Paths:\n"); \
    .libPaths() |> cat(sep="\n"); \
    if (!require("kinfitr")) stop("Failed to install kinfitr package")'

# Expose port for RStudio Server
EXPOSE 8787

# Start RStudio Server
CMD ["/init"]