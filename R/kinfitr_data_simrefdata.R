#' Simulated Reference Tissue Dataset
#'
#' A dataset containing simulated time activity curves (TACs) and blood
#' data for reference tissue models.
#'
#' @format A nested data frame with 24 rows and 4 variables:
#' \describe{
#'   \item{Subjname}{Descriptor of which individual was measured}
#'   \item{PETNo}{Descriptor of the PET number}
#'   \item{PET}{Descriptor of the PET measurement}
#'   \item{tacs}{Nested tibble of the TAC data. Consists of Times, in minutes; Weights, between 0 and 1; StartTime, in minutes; Duration, in minutes, and TACs for a simulated reference region, and three target regions.}
#' }
"simref"
#' @usage data("sim_ref")
#' @keywords PET, Reference, TAC
#' @docType data
