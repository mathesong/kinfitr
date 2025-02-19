# set url_base to http://openneuropet.org for deployment
# url_base <- "http://54.144.240.214"
# get_url <- paste0(url_base, "/check/kinfitr/")
# post_url <- paste0(url_base, "/kinfitr/")



#' Get Telemetry Data
#'
#' Checks to see what has been posted to the URL endpoint. Should return
#' location and any any other data that gets put there with send_telemetry
#'
#' @param url The url endpoint
#' @param number_of_records How many records. 0 for all.
#'
#' @returns A list of the recorded telemetry data
get_telemetry <- function(url = "http://54.144.240.214/check/kinfitr/",
                          number_of_records = 0) {


  req <- httr2::request(paste0(url, as.character(number_of_records))) %>%
    httr2::req_headers("Accept" = "application/json") %>%
    httr2::req_retry(max_tries = 3) %>%
    httr2::req_timeout(3)

  response <- httr2::req_perform(req)

  httr2::resp_body_json(response)

}

#' Sending telemetry data
#'
#' Sending telemetry data to the URL endpoint
#'
#' @param telemetry_json_data The JSON data content to be sent
#' @param url The url endpoint
#'
#' @returns The status code
send_telemetry <- function(telemetry_json_data,
                           url = "http://54.144.240.214/kinfitr/") {

  no_track <- Sys.getenv("KINFITR_NO_TRACK")

  if (tolower(no_track) == "true") {
    return(NULL)
  } else {

    try(
      {
        req <- httr2::request(url) %>%
          httr2::req_retry(max_tries = 3) %>%
          httr2::req_timeout(3)

        response <- req %>%
          httr2::req_body_json(data = telemetry_json_data) %>%
          httr2::req_perform()

        return_values <- list(status_code = response$status_code)

        return(return_values)
      },
      silent = FALSE
    )
  }
}
