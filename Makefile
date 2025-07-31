# Variables
RSTUDIO_IMAGE = kinfitr-rstudio.sif
DEV_IMAGE = kinfitr-dev.sif
BIND_PATH = ${PWD}
DOCKER_DEV_IMAGE = kinfitr-dev-docker
DOCKER_RSTUDIO_IMAGE = kinfitr-rstudio-docker

# Default target
all: build-all

# Build all containers
build-all: build-rstudio build-dev build-docker-dev build-docker-rstudio

# Build RStudio container
build-rstudio:
	apptainer build containers/$(RSTUDIO_IMAGE) containers/rstudio.Singularity

# Build Development container
build-dev:
	apptainer build containers/$(DEV_IMAGE) containers/dev.Singularity

# Build Docker development container for VSCode
build-docker-dev:
	docker build -t containers/$(DOCKER_DEV_IMAGE) -f Dockerfile.dev .

# Build Docker RStudio container for VSCode
build-docker-rstudio:
	docker build -t containers/$(DOCKER_RSTUDIO_IMAGE) -f containers/rstudio.Dockerfile .

# Run RStudio container
run-rstudio:
	apptainer run --bind $(BIND_PATH):/home/rstudio/kinfitr containers/$(RSTUDIO_IMAGE)

# Run Development container
run-dev:
	apptainer shell --bind $(BIND_PATH):/workspace containers/$(DEV_IMAGE)

# Run Docker development container
run-docker-dev:
	docker run -it --rm \
		-v $(BIND_PATH):/workspace \
		$(DOCKER_DEV_IMAGE)

# Clean up built images
clean:
	rm -f containers/$(RSTUDIO_IMAGE) containers/$(DEV_IMAGE)
	docker rmi $(DOCKER_DEV_IMAGE) || true

.PHONY: all build-all build-rstudio build-dev build-docker-dev run-rstudio run-dev run-docker-dev clean 