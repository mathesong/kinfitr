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

write_installation_info <- function(install_info,
                                    path_to_installation_info_json) {
  # Create a list to write to JSON
  install_info$installed_on <- Sys.Date()

  # Write the list to the JSON file
  tryCatch(
    {
      jsonlite::write_json(install_info, path_to_installation_info_json)
    },
    error = function(e) {
      # Fail silently
      message(
        "Failed to write installation info to: ", path_to_installation_info_json
      )
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
  path_to_installation_info_json <- file.path(
    libname, "kinfitr", "installation_settings.json"
  )
  install_info <- read_installation_info(
    path_to_installation_info_json
  )
  # Skip telemetry during tests
  if (is_loading_for_tests()) {
    return(NULL)
  }

  if (isTRUE(getOption("kinfitr.no_track") |> as.logical()) ||
      isTRUE(Sys.getenv("KINFITR_NO_TRACK") |> as.logical()) ||
      isTRUE(install_info$no_track |> as.logical())) {

    install_info$no_track <- TRUE
    write_installation_info(install_info, path_to_installation_info_json)

    return(NULL)
  }

  if (!isTRUE(getOption("kinfitr.no_track") |> as.logical()) &&
      !isTRUE(Sys.getenv("KINFITR_NO_TRACK") |> as.logical()) &&
      !isTRUE(install_info$no_track |> as.logical()) &&
      !isTRUE(is_loading_for_tests())) {

    send_telemetry(list("kinfitr_usage" = "package_installed"))
    write_installation_info(install_info, path_to_installation_info_json)
  }

  # display this message at install
  if (!isTRUE(install_info$displayed_message %>% as.logical())) {
    message(
      paste(
        "Please note that kinfitr will send telemetry information to the",
        "kinfitr developers in order to track real world usage, which is",
        "critical for obtaining funding to continue future development.",
        "This consists of information about approximate location and operating",
        "system.\n\n",
        "To opt-out of sending this information and disable usage tracking, set",
        "options(kinfitr.no_track = TRUE) or the environment variable",
        "KINFITR_NO_TRACK to TRUE. To hide this message, set",
        "options(kinfitr.no_track = FALSE) or the environment variable",
        "KINFITR_NO_TRACK to FALSE .",
        sep = " "
      )
    )
    install_info$displayed_message <- TRUE
    write_installation_info(install_info, path_to_installation_info_json)
  }
}

.onLoad <- function(libname, pkgname) {
  path_to_installation_info_json <- file.path(
    libname, "kinfitr", "installation_settings.json"
  )

  install_info <- read_installation_info(
    path_to_installation_info_json
  )
  # check to see if the user has opted out of tracking at a
  # system environment, R environment, or via config file
  # or if this library is being tested
  if (isTRUE(getOption("kinfitr.no_track") |> as.logical()) ||
      isTRUE(Sys.getenv("KINFITR_NO_TRACK") |> as.logical()) ||
      isTRUE(install_info$no_track |> as.logical())) {

    install_info$no_track <- TRUE
    write_installation_info(install_info, path_to_installation_info_json)

    return(NULL)
  }

  if (isTRUE(is_loading_for_tests())) {
    # don't do anything if running tests
    return(NULL)
  }

  if (!isTRUE(getOption("kinfitr.no_track") |> as.logical()) &&
      !isTRUE(Sys.getenv("KINFITR_NO_TRACK") |> as.logical()) &&
      !isTRUE(install_info$no_track |> as.logical()) &&
      !isTRUE(is_loading_for_tests())) {

    # Send telemetry when package is loaded
    send_telemetry(list("kinfitr_usage" = "package_loaded"))

    # if none of the conditionals are met in the first statement, then
    # set tracking to true (no_track=FALSE)
    install_info$no_track <- FALSE
    write_installation_info(install_info, path_to_installation_info_json)
  }
}
