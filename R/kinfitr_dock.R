# Cache environment for DOCK database queries
.dock_cache <- new.env(parent = emptyenv())


# Internal: download and extract DOCK database from GitHub
# Downloads the repo tarball once per session and caches locally
.dock_download <- function() {

  if (exists("dock_dir", envir = .dock_cache)) {
    dock_dir <- get("dock_dir", envir = .dock_cache)
    if (dir.exists(dock_dir)) return(dock_dir)
  }

  url <- paste0(
    "https://github.com/",
    "Database-of-CNS-PET-Kinetic-parameters/DOCK-PET/",
    "archive/refs/heads/main.tar.gz"
  )

  tmp_tar <- tempfile(fileext = ".tar.gz")
  tmp_dir <- tempfile(pattern = "dock_")

  tryCatch({
    utils::download.file(url, tmp_tar, quiet = TRUE, mode = "wb")
    dir.create(tmp_dir)
    utils::untar(tmp_tar, exdir = tmp_dir)
    unlink(tmp_tar)
  },
    error = function(e) {
      unlink(tmp_tar)
      unlink(tmp_dir, recursive = TRUE)
      stop("Failed to download DOCK database from GitHub. ",
           "Check your internet connection.\n",
           "Original error: ", e$message,
           call. = FALSE)
    }
  )

  dock_dir <- file.path(
    tmp_dir, "DOCK-PET-main", "CNS-PET tracer JSON database"
  )

  if (!dir.exists(dock_dir)) {
    unlink(tmp_dir, recursive = TRUE)
    stop("Unexpected DOCK repository structure.", call. = FALSE)
  }

  assign("dock_dir", dock_dir, envir = .dock_cache)
  dock_dir
}


# Internal: get list of all tracers from DOCK database
# Returns a tibble with columns: tracer, target, path
.dock_tracer_list <- function() {

  if (exists("tracer_list", envir = .dock_cache)) {
    return(get("tracer_list", envir = .dock_cache))
  }

  dock_dir <- .dock_download()

  json_files <- list.files(dock_dir, pattern = "\\.json$",
                           recursive = TRUE, full.names = TRUE)

  tracer_list <- tibble::tibble(
    tracer = tools::file_path_sans_ext(basename(json_files)),
    target = basename(dirname(json_files)),
    path   = json_files
  )

  assign("tracer_list", tracer_list, envir = .dock_cache)
  tracer_list
}


# Internal: read and parse a single tracer JSON from local cache
.dock_fetch_tracer <- function(path) {
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}


# Internal: recursively extract region-level mean/SD from nested parameter list
.dock_extract_regions <- function(x, param_name) {

  # Check if this level has _mean and _SD keys (leaf node)
  keys <- names(x)
  mean_key <- grep("_mean$", keys, value = TRUE)
  sd_key   <- grep("_SD$", keys, value = TRUE)

  if (length(mean_key) == 1 && length(sd_key) == 1) {
    region <- sub("_mean$", "", mean_key)
    if (nchar(region) == 0) return(NULL)  # skip empty region names
    return(tibble::tibble(
      parameter = param_name,
      region    = region,
      mean      = if (is.null(x[[mean_key]])) NA_real_ else as.numeric(x[[mean_key]]),
      sd        = if (is.null(x[[sd_key]]))   NA_real_ else as.numeric(x[[sd_key]])
    ))
  }

  # Otherwise recurse into sub-elements that are lists
  results <- list()
  for (nm in keys) {
    if (is.list(x[[nm]])) {
      results <- c(results, list(.dock_extract_regions(x[[nm]], param_name)))
    }
  }

  if (length(results) > 0) {
    dplyr::bind_rows(results)
  } else {
    NULL
  }
}


# Internal: parse full JSON into a flat tibble
.dock_parse_parameters <- function(json, tracer_name) {

  # Handle key name inconsistencies
  drug_key <- intersect(
    c("Drug data", "Drug information"),
    names(json)
  )[1]
  journal_key <- intersect(
    c("Journal data", "Journal Information", "Journal information"),
    names(json)
  )[1]
  study_key <- intersect(
    c("Study data", "Study information"),
    names(json)
  )[1]
  param_key <- intersect(
    c("Parameters", "Parameters ", "Parameters information"),
    names(json)
  )[1]

  # Extract metadata
  drug <- json[[drug_key]][[1]]
  target <- if (!is.null(drug$target)) drug$target else NA_character_

  journal <- json[[journal_key]][[1]]
  author  <- if (!is.null(journal$auther)) journal$auther else NA_character_
  year    <- if (!is.null(journal$year))   journal$year   else NA_character_
  jrnl    <- if (!is.null(journal$journal)) journal$journal else NA_character_
  doi     <- if (!is.null(journal$DOI))    journal$DOI    else NA_character_

  study <- json[[study_key]][[1]]
  n_subjects <- if (!is.null(study[["no.subjects"]])) {
    as.integer(study[["no.subjects"]])
  } else {
    NA_integer_
  }

  # Extract parameters

  params <- json[[param_key]]
  if (is.null(params) || length(params) < 2) {
    warning("No parameter data found for tracer: ", tracer_name, call. = FALSE)
    return(tibble::tibble())
  }

  # First element is method info
  method_info <- params[[1]]
  method <- if (!is.null(method_info$method)) method_info$method else NA_character_
  ref_region <- if (!is.null(method_info[["reference region"]])) {
    method_info[["reference region"]]
  } else {
    NA_character_
  }

  # Second element contains parameter data
  param_data <- params[[2]]
  param_names <- names(param_data)

  all_regions <- list()
  for (pname in param_names) {
    if (is.list(param_data[[pname]])) {
      extracted <- .dock_extract_regions(param_data[[pname]], pname)
      if (!is.null(extracted) && nrow(extracted) > 0) {
        all_regions <- c(all_regions, list(extracted))
      }
    }
  }

  if (length(all_regions) == 0) {
    warning("No parameter values found for tracer: ", tracer_name, call. = FALSE)
    return(tibble::tibble())
  }

  result <- dplyr::bind_rows(all_regions)

  # Add metadata columns
  result$tracer           <- tracer_name
  result$target           <- target
  result$method           <- method
  result$reference_region <- ref_region
  result$author           <- author
  result$year             <- year
  result$journal          <- jrnl
  result$doi              <- doi
  result$n_subjects       <- n_subjects

  # Reorder columns
  result <- result[, c(
    "tracer", "target", "method", "reference_region",
    "parameter", "region", "mean", "sd",
    "author", "year", "journal", "doi", "n_subjects"
  )]

  result
}


#' Find tracers in the DOCK database
#'
#' Searches the DOCK (Database of CNS PET Kinetic parameters) database for
#' tracers matching a search string. Requires an internet connection to query
#' the database from GitHub.
#'
#' @param tracer_grob A character string to search for among tracer names
#'   (case-insensitive). For example, \code{"raclopride"}, \code{"PBR"},
#'   or \code{"WAY"}.
#'
#' @return A tibble with columns \code{tracer} and \code{target}, containing
#'   all matching tracers. Returned invisibly.
#'
#' @details The DOCK database is described in Kudo et al. (2024). Database of
#'   CNS PET kinetic parameters (DOCK): a new database of average
#'   parameter estimates for neuroreceptors and transporters. Annals of
#'   Nuclear Medicine. DOI: 10.1007/s12149-024-01947-z
#'
#' @examples
#' \dontrun{
#' dock_find_tracer("raclopride")
#' dock_find_tracer("PBR")
#' dock_find_tracer("WAY")
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
dock_find_tracer <- function(tracer_grob) {

  tracer_list <- .dock_tracer_list()

  matches <- grep(tracer_grob, tracer_list$tracer, ignore.case = TRUE)

  if (length(matches) == 0) {
    # Try approximate matching
    matches <- agrep(tracer_grob, tracer_list$tracer, ignore.case = TRUE)
    if (length(matches) > 0) {
      message("No exact matches found. Approximate matches:")
    } else {
      message("No tracers found matching '", tracer_grob, "'")
      return(invisible(tibble::tibble(tracer = character(), target = character())))
    }
  }

  result <- tracer_list[matches, c("tracer", "target")]

  for (i in seq_len(nrow(result))) {
    message(result$tracer[i], "  [", result$target[i], "]")
  }

  invisible(result)
}


#' Query kinetic parameters from the DOCK database
#'
#' Queries the DOCK (Database of CNS PET Kinetic parameters) database for
#' kinetic parameter estimates for a given radioligand. Requires an internet
#' connection to fetch data from GitHub.
#'
#' @param tracer_name The name of the tracer to query, e.g.
#'   \code{"11C_raclopride"} or \code{"18F_FDG"}. Use
#'   \code{\link{dock_find_tracer}} to search for available tracer names.
#' @param model Optional character string to filter by kinetic model
#'   (case-insensitive partial match), e.g. \code{"2TCM"}.
#' @param region Optional character string to filter by brain region
#'   (case-insensitive partial match), e.g. \code{"Thalamus"}.
#' @param parameter Optional character string to filter by parameter name
#'   (case-insensitive partial match), e.g. \code{"Vt"}, \code{"BPND"},
#'   \code{"K1"}.
#'
#' @return A tibble with columns: \code{tracer}, \code{target}, \code{method},
#'   \code{reference_region}, \code{parameter}, \code{region}, \code{mean},
#'   \code{sd}, \code{author}, \code{year}, \code{journal}, \code{doi},
#'   \code{n_subjects}. Rows where both \code{mean} and \code{sd} are
#'   \code{NA} are excluded.
#'
#' @details The DOCK database is described in Kudo et al. (2024). Database of
#'   CNS PET kinetic parameters (DOCK): a new database of average
#'   parameter estimates for neuroreceptors and transporters. Annals of
#'   Nuclear Medicine. DOI: 10.1007/s12149-024-01947-z
#'
#' @examples
#' \dontrun{
#' dock_query_tracer("11C_raclopride")
#' dock_query_tracer("11C_DASB", parameter = "Vt")
#' dock_query_tracer("11C_DASB", region = "Thalamus", parameter = "K1")
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export
dock_query_tracer <- function(tracer_name, model = NULL,
                              region = NULL, parameter = NULL) {

  tracer_list <- .dock_tracer_list()

  # Find exact match

  idx <- which(tracer_list$tracer == tracer_name)

  if (length(idx) == 0) {
    # Try case-insensitive exact match
    idx <- which(tolower(tracer_list$tracer) == tolower(tracer_name))
  }

  if (length(idx) == 0) {
    # Try partial match and suggest
    partial <- grep(tracer_name, tracer_list$tracer, ignore.case = TRUE)
    if (length(partial) > 0) {
      suggestions <- paste(tracer_list$tracer[partial], collapse = ", ")
      stop("Tracer '", tracer_name, "' not found. Did you mean: ",
           suggestions, "?",
           "\nUse dock_find_tracer() to search for available tracers.",
           call. = FALSE)
    } else {
      stop("Tracer '", tracer_name, "' not found in the DOCK database. ",
           "Use dock_find_tracer() to search for available tracers.",
           call. = FALSE)
    }
  }

  # Fetch and parse
  path <- tracer_list$path[idx[1]]
  json <- .dock_fetch_tracer(path)
  result <- .dock_parse_parameters(json, tracer_list$tracer[idx[1]])

  if (nrow(result) == 0) return(result)

  # Remove rows where both mean and sd are NA

  result <- result[!(is.na(result$mean) & is.na(result$sd)), ]

  # Apply optional filters
  if (!is.null(model)) {
    result <- result[grep(model, result$method, ignore.case = TRUE), ]
  }
  if (!is.null(region)) {
    result <- result[grep(region, result$region, ignore.case = TRUE), ]
  }
  if (!is.null(parameter)) {
    result <- result[grep(parameter, result$parameter, ignore.case = TRUE), ]
  }

  result
}
