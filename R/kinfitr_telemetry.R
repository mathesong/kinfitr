
get_url <- "http://54.144.240.214/check/kinfitr/"
post_url <- "http://54.144.240.214/kinfitr/"

get_telemetry <- function(url = get_url, number_of_records = 0) {
  # checks to see what's been posted to the url endpoint, should return location,
  # and any other data that gets put there with send_telemetry
  req <- httr2::request(paste(url, as.character(number_of_records), sep = ""))
  req |> httr2::req_headers("Accept" = "application/json")
  req |> httr2::req_retry(max_tries = 3)
  req |> httr2::req_timeout(3)
  response <- httr2::req_perform(req)
  return(resp_body_json(response))
}

#' Send telemetry data
#' @keywords internal
send_telemetry <- function(telemetry_json_data, url = post_url) {
  # Add debug messages
  message("Debug: Starting send_telemetry")
  message("Debug: Telemetry data: ", toString(telemetry_json_data))

  no_track <- Sys.getenv("KINFITR_NO_TRACK")
  message("Debug: KINFITR_NO_TRACK value: ", no_track)

  if (tolower(no_track) == "true") {
    message("Debug: Tracking disabled")
    return(NULL)
  } else {
    message("Debug: Attempting to send telemetry")
    try(
      {
        req <- httr2::request(url)
        # Fix the request chain - these need to be assigned
        req <- httr2::req_retry(req, max_tries = 3)
        req <- httr2::req_timeout(req, 3)

        message("Debug: Sending to URL: ", url)
        response <- req |>
          httr2::req_body_json(data = telemetry_json_data) |>
          httr2::req_perform()

        message("Debug: Response status: ", response$status_code)
        return_values <- list(status_code = response$status_code)
        return(return_values)
      },
      silent = FALSE # Changed to FALSE to see any errors
    )
  }
}
