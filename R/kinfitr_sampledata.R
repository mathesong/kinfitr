#' PBR28 Test-Retest Data
#'
#' A dataset containing the time activity curves (TACs) and blood
#' data from a study of [11C]PBR28.
#'
#' @format A nested data frame with 24 rows and 10 variables:
#' \describe{
#'   \item{PET}{Descriptor of the PET measurement}
#'   \item{Subjname}{Descriptor of which individual was measured}
#'   \item{PETNo}{Descriptor of the PET number}
#'   \item{tacs}{Nested tibble of the TAC data. Consists of Times, in seconds; Weights, between 0 and 1; StartTime, in seconds; Duration, in seconds, and TACs for Frontal Cortex (FC), Temporal Cortex (TC), Striatum (STR), Thalamus (THA), Whole Brain (WB), and Cerebellum (CBL)}
#'   \item{procblood}{Nested tibble of the processed blood data which was generated for Matheson et al. (2017). Preprocessing was performed using in-house methods. This data does not contain information about automatic and manual sampling. Consists of Time in seconds (Time), blood radioactivity following dispersion correction (Cbl_dispcorr) and plasma radioactivity following metabolite correction (Cpl_metabcorr)}
#'   \item{Genotype}{Refers to the High Affinity Binder (HAB) or Mixed Affinity Binder (MAB) genotype}
#'   \item{injRad}{Injected radioactivity}
#'   \item{jsondata}{This information is the contents of the BIDS json sidecar for each PET measurement. This contains all the raw blood data.}
#' }
#' @source \url{https://ejnmmires.springeropen.com/articles/10.1186/s13550-017-0304-1}
#' @references Matheson, Granville J, Plav\'en-Sigray, Pontus, Forsberg, Anton, Varrone, Andrea, Farde, Lars & Cervenka, Simon (2017). Assessment of simplified ratio-based approaches for quantification of PET [11 C] PBR28 data. EJNMMI research, 7, 58.
"pbr28"
#' @usage data("pbr28")
#' @keywords PET, TSPO, PBR28, TAC
#' @docType data
