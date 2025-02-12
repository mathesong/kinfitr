# Variables
RSTUDIO_IMAGE = kinfitr-rstudio.sif
DEV_IMAGE = kinfitr-dev.sif
BIND_PATH = ${PWD}

# Default target
all: build-all

# Build all containers
build-all: build-rstudio build-dev

# Build RStudio container
build-rstudio:
	apptainer build $(RSTUDIO_IMAGE) containers/rstudio.Singularity

# Build Development container
build-dev:
	apptainer build $(DEV_IMAGE) containers/dev.Singularity

# Run RStudio container
run-rstudio:
	apptainer run --bind $(BIND_PATH):/home/rstudio/kinfitr $(RSTUDIO_IMAGE)

# Run Development container
run-dev:
	apptainer shell --bind $(BIND_PATH):/workspace $(DEV_IMAGE)

# Clean up built images
clean:
	rm -f $(RSTUDIO_IMAGE) $(DEV_IMAGE)

.PHONY: all build-all build-rstudio build-dev run-rstudio run-dev clean 