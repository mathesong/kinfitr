#!/bin/bash

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to clean up on exit
cleanup() {
    echo -e "\n${YELLOW}Cleaning up...${NC}"
    docker stop mathesong-kinfitr 2>/dev/null
    exit 0
}

# Set up trap for cleanup on script exit
trap cleanup SIGINT SIGTERM

# Function to check if Docker is running
check_docker() {
    if ! docker info >/dev/null 2>&1; then
        echo -e "${RED}Error: Docker is not running${NC}"
        exit 1
    fi
}

# Function to check if port 8787 is available
check_port() {
    if lsof -Pi :8787 -sTCP:LISTEN -t >/dev/null ; then
        echo -e "${RED}Error: Port 8787 is already in use${NC}"
        exit 1
    fi
}

# Function to check if RStudio is ready
wait_for_rstudio() {
    echo -e "${YELLOW}Waiting for RStudio to be ready...${NC}"
    local max_attempts=30
    local attempt=1
    
    while ! curl -s http://localhost:8787 > /dev/null; do
        if [ $attempt -eq $max_attempts ]; then
            echo -e "${RED}Error: RStudio failed to start${NC}"
            cleanup
            exit 1
        fi
        echo -n "."
        sleep 1
        ((attempt++))
    done
    echo -e "\n${GREEN}RStudio is ready!${NC}"
}

# Main script
echo -e "${YELLOW}Starting kinfitr RStudio environment...${NC}"

# Check prerequisites
check_docker
check_port

# Build the Docker image if it doesn't exist
if [[ "$(docker images -q mathesong/kinfitr:latest 2> /dev/null)" == "" ]]; then
    echo -e "${YELLOW}Building Docker image...${NC}"
    if ! docker build -t mathesong/kinfitr -f Dockerfile-Rstudio .; then
        echo -e "${RED}Error: Failed to build Docker image${NC}"
        exit 1
    fi
fi

# Launch the Docker container in the background
echo -e "${YELLOW}Starting RStudio container...${NC}"
if ! docker run --rm -d -p 8787:8787 --name mathesong-kinfitr mathesong/kinfitr:latest; then
    echo -e "${RED}Error: Failed to start container${NC}"
    exit 1
fi

# Wait for RStudio to be ready
wait_for_rstudio

# Open browser based on OS
echo -e "${YELLOW}Opening RStudio in browser...${NC}"
case "$(uname -s)" in
    Linux*)     xdg-open http://localhost:8787 ;;
    Darwin*)    open http://localhost:8787 ;; # macOS
    MINGW*|CYGWIN*|MSYS*)    start http://localhost:8787 ;; # Windows
    *)          echo -e "${YELLOW}Please open http://localhost:8787 in your browser${NC}" ;;
esac

echo -e "${GREEN}
RStudio is now running at http://localhost:8787
Press Ctrl+C to stop the container and cleanup
${NC}"

# Wait for Ctrl+C
docker logs -f mathesong-kinfitr 