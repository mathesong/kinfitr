library(httr2)

get_url <- "http://127.0.0.1:8000/check/kinfitr/"
post_url <- "http://127.0.0.1:8000/kinfitr/"

get_telemetry <- function(url = get_url, number_of_records = 0) {
  # checks to see what's been posted to the url endpoint, should return location,
  # and any other data that gets put there with send_telemetry
  req <- request(paste(url, as.character(number_of_records), sep = ""))
  req |> req_headers("Accept" = "application/json")
  req |> req_retry(max_tries = 3)
  req |> req_timeout(3)
  response <- req_perform(req)
  return(resp_body_json(response))
}

#' Send telemetry data
#' @keywords internal
send_telemetry <- function(telemetry_json_data, url = post_url) {
  # Accepts a list of telemetry data (worth transmitting of course)
  # and converts it such that it can be sent as a json to the endpoint.
  # Additionally, checks to to see if environment variable
  # KINFITR_NO_TRACK is set if it is, this does nothing.

  no_track <- Sys.getenv("KINFITR_NO_TRACK")

  if (tolower(no_track) == "true") {
    return(NULL)
  } else {
    try(
      {
        req <- request(url)
        req |> req_retry(max_tries = 3)
        req |> req_timeout(3)
        response <- req |>
          req_body_json(data = telemetry_json_data) |>
          req_perform()
        return_values <- list(status_code = response$status_code)
        return(return_values)
      },
      silent = TRUE
    )
  }
}


# list_to_send <- list("part1" = "one", "part2" = "new")
# end_telemetry(list_to_send)
