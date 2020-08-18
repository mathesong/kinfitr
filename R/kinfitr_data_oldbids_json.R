#' Old BIDS specification JSON Data
#'
#' A dataset containing the old way that data was structured, for managing
#' the transition to the new PET BIDS structure. This is from the PBR28 data.
#'
#' @format A nested list of the way that data used to be specified in the PET BIDS format:
#' \describe{
#'   \item{Info}{PET measurement information}
#'   \item{Time}{Time information}
#'   \item{Radiochem}{Radiochemistry information}
#'   \item{Plasma}{Plasma data}
#'   \item{Metabolite}{Metabolite data}
#'   \item{Blood}{Blood data}
#' }
#' @source \url{https://ejnmmires.springeropen.com/articles/10.1186/s13550-017-0304-1}
#' @references Matheson, Granville J, Plav\'en-Sigray, Pontus, Forsberg, Anton, Varrone, Andrea, Farde, Lars & Cervenka, Simon (2017). Assessment of simplified ratio-based approaches for quantification of PET [11 C] PBR28 data. EJNMMI research, 7, 58.
"pbr28"
#' @usage data("oldbids_json")
#' @keywords PET, TSPO, PBR28, JSON
#' @docType data
