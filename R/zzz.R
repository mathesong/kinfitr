.onAttach <- function(libname, pkgname) {
  message("Loading kinfitr")
  message("new debug message")
  # Display opt-out message
  if (
      !isTRUE(getOption("kinfitr.no_track")) &&
        Sys.getenv("KINFITR_NO_TRACK") != "TRUE") {
    message(
      "Opt-out of sending tracking information to the KinFitR developers. ",
      "This information provides an indicator of real world usage crucial for ",
      "obtaining funding. To disable tracking set ",
      "options(kinfitr.no_track = TRUE) or the environment variable ",
      "KINFITR_NO_TRACK to TRUE), to hide this message set ",
      "options(kinfitr.no_track = FALSE) or the environment variable ",
      "KINFITR_NO_TRACK to FALSE"
    )
  }
  # # Check if user has opted out
  if (
      isTRUE(getOption("kinfitr.no_track")) ||
        isTRUE(Sys.getenv("KINFITR_NO_TRACK"))) {
    # do nothing
    NULL
  } else {
    # Send telemetry when package is loaded
    message("Sending telemetry")
    send_telemetry(list("kinfitr_usage" = "package_loaded"))
  }
}
