library(httr2)

get_url <- "http://127.0.0.1:8000/check/kinfitr/"
post_url <- "http://127.0.0.1:8000/kinfitr/"

get_telemetry <- function(url = get_url, number_of_records = 0) {
  # checks to see what's been posted to the url endpoint, should return
  # location and any other data that gets put there with send_telemetry
  req <- request(paste(url, as.character(number_of_records), sep = ""))
  req |> req_headers("Accept" = "application/json")
  req |> req_retry(max_tries = 3)
  req |> req_timeout(3)
  response <- req_perform(req)
  return(resp_body_json(response))
}

.onLoad <- function(libname, pkgname) {
  # Display opt-out message
  if (
      isFALSE(getOption("kinfitr.no_track")) ||
        isFALSE(Sys.getenv("KINFITR_NO_TRACK"))) {
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
  # Check if user has opted out
  if (
      isTRUE(getOption("kinfitr.no_track")) ||
        isTRUE(Sys.getenv("KINFITR_NO_TRACK"))) {
    # Send telemetry when package is loaded
    send_telemetry(list("kinfitr_usage" = "package_loaded"))
  }
}
