

read_installation_info <- function(path_to_installation_info_json) {
  # Initialize a list to store the key-value pairs
  installation_info <- list(
    package = NULL,
    installed_on = NULL,
    no_track = NULL,
    displayed_message = NULL
  )

  # Read the JSON file
  tryCatch(
    {
      info <- jsonlite::fromJSON(path_to_installation_info_json)
      # update installation info with what's in the json file
      installation_info <- as.list(info)
    },
    error = function(e) {
      # Fail silently, keep default NULL values
    }
  )

  return(installation_info)
}

write_installation_info <- function(
    install_info, path_to_installation_info_json) {
  # Create a list to write to JSON
  install_info$installed_on <- Sys.Date()

  # Write the list to the JSON file
  tryCatch(
    {
      jsonlite::write_json(install_info, path_to_installation_info_json)
      message("Installation settings written to: ", path_to_installation_info_json)
    },
    error = function(e) {
      # Fail silently
      message("Failed to write installation info to: ", path_to_installation_info_json)
      message("Error: ", e)
    }
  )
}

is_loading_for_tests <- function() {
  !interactive() && (
    identical(Sys.getenv("DEVTOOLS_LOAD"), "kinfitr") ||
      identical(Sys.getenv("TESTTHAT"), "true")
  )
}

.onAttach <- function(libname, pkgname) {
  path_to_installation_info_json <- file.path(libname, "kinfitr", "installation_settings.json")
  # Check if tracking is disabled
  if (!isTRUE(getOption("kinfitr_no_track"))) {
    message("Opt-out of sending tracking information to
    the KinFitR developers.")
  }

  # Skip telemetry during tests
  if (is_loading_for_tests()) {
    return(NULL)
  }

  install_info <- read_installation_info(path_to_installation_info_json)
  if (!isTRUE(getOption("kinfitr.no_track")) &&
        Sys.getenv("KINFITR_NO_TRACK") != "TRUE" &&
        !isTRUE(install_info$displayed_message)) {
    message(
      "Opt-out of sending tracking information to the KinFitR developers. ",
      "This information provides an indicator of real world usage crucial for ",
      "obtaining funding. To disable tracking set ",
      "options(kinfitr.no_track = TRUE) or the environment variable ",
      "KINFITR_NO_TRACK to TRUE), to hide this message set ",
      "options(kinfitr.no_track = FALSE) or the environment variable ",
      "KINFITR_NO_TRACK to FALSE"
    )
    send_telemetry(list("kinfitr_usage" = "package_installed"))
    install_info$displayed_message <- TRUE
    message("Writing installation info to: ", path_to_installation_info_json)
    write_installation_info(install_info, path_to_installation_info_json)
  }
}

.onLoad <- function(libname, pkgname) {
  path_to_installation_info_json <- file.path(libname, "kinfitr", "installation_settings.json")
  if (isTRUE(getOption("kinfitr.no_track")) ||
        Sys.getenv("KINFITR_NO_TRACK") == "TRUE" ||
        isTRUE(is_loading_for_tests())) {
    # do nothing
    install_info <- read_installation_info(path_to_installation_info_json)
    install_info$no_track <- TRUE
    write_installation_info(install_info, path_to_installation_info_json)
    message("Installation info written to: ", path_to_installation_info_json)
    return(NULL)
  } else {
    # Send telemetry when package is loaded
    send_telemetry(list("kinfitr_usage" = "package_loaded"))
  }

}
