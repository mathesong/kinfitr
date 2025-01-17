# Set working directory based on environment
if (Sys.getenv("DOCKER") == "true") {
    # Docker environment - use absolute path
    setwd("/home/rstudio/kinfitr")
    # Add package to search path
    .First <- function() {
        message("Loading kinfitr from source...")
        devtools::load_all(".")
        message("✓ Loaded kinfitr")
    }
} else {
    # Local environment - use relative path from project root
    # Get the project root directory (where .Rprofile lives)
    project_root <- normalizePath(dirname(sys.frame(1)$ofile))
    setwd(project_root)
    # Add package to search path
    .First <- function() {
        message("Loading kinfitr from source...")
        devtools::load_all(".")
        message("✓ Loaded kinfitr")
    }
}

# Auto-load packages from DESCRIPTION file when in RStudio
if (interactive() && Sys.getenv("RSTUDIO") == "1") {
    message("Loading project dependencies from DESCRIPTION file...")

    # Function to safely load packages from DESCRIPTION
    load_description_packages <- function() {
        if (file.exists("DESCRIPTION")) {
            deps <- read.dcf("DESCRIPTION", fields = c("Imports", "Depends"))
            packages <- unlist(strsplit(paste(deps[!is.na(deps)], collapse = ","), ","))
            packages <- gsub("\\s+", "", packages)
            packages <- packages[packages != "R"]

            # Load each package
            for (package in packages) {
                # Remove version specifications if present
                package <- gsub("\\s*\\(.*\\)", "", package)
                tryCatch({
                    library(package, character.only = TRUE)
                    message(sprintf("✓ Loaded %s", package))
                }, error = function(e) {
                    warning(sprintf("Failed to load %s: %s", package, e$message))
                })
            }
        } else {
            warning("No DESCRIPTION file found in current directory")
        }
    }

    load_description_packages()
}
