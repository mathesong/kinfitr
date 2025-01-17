# Set working directory based on environment
if (Sys.getenv("DOCKER") == "true") {
    # Docker environment - use absolute path
    setwd("/home/rstudio/kinfitr")
    # Add R package directory to the search path
    .libPaths(c("/home/rstudio/kinfitr/R", .libPaths()))
} else {
    # Local environment - use relative path from project root
    project_root <- normalizePath(dirname(sys.frame(1)$ofile))
    setwd(project_root)
    .libPaths(c(file.path(project_root, "R"), .libPaths()))
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

    # Source all R files directly
    message("Loading kinfitr source files...")
    r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
    for (file in r_files) {
        tryCatch({
            source(file)
            message(sprintf("✓ Loaded %s", basename(file)))
        }, error = function(e) {
            warning(sprintf("Failed to load %s: %s", basename(file), e$message))
        })
    }
}
