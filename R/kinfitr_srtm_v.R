#' Simplified Reference Tissue Model with Blood Volumes
#'
#' Function to fit the SRTM_V model of Tomasi et al (2008) to data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the
#'   time halfway through the frame as well as a zero. If a time zero frame is
#'   not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference
#'   tissue for each frame. We include zero at time zero: if not included, it is
#'   added.
#' @param roitac Numeric vector of radioactivity concentrations in the target
#'   tissue for each frame. We include zero at time zero: if not included, it is
#'   added.
#' @param bloodtac Numeric vector of radioactivity concentrations in the blood
#'   for each frame. We include zero at time zero: if not included, it is added.
#' @param weights Optional. Numeric vector of the weights assigned to each frame
#'   in the fitting. We include zero at time zero: if not included, it is added.
#'   If not specified, uniform weights will be used.
#' @param vBr Optional. The blood volume fraction of the reference region.  If
#'   not specified, this will be fitted. This parameter was fixed in the
#'   original article.
#' @param frameStartEnd Optional. This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). This is to assess time
#'   stability.
#' @param R1.start Optional. Starting parameter for fitting of R1. Default is 1.
#' @param R1.lower Optional. Lower bound for the fitting of R1. Default is
#'   0.0001.
#' @param R1.upper Optional. Upper bound for the fitting of R1. Default is 10.
#' @param k2.start Optional. Starting parameter for fitting of k2. Default is
#'   0.1.
#' @param k2.lower Optional. Lower bound for the fitting of k2. Default is
#'   0.0001.
#' @param k2.upper Optional. Upper bound for the fitting of k2. Default is 1.
#' @param bp.start Optional. Starting parameter for fitting of bp. Default is
#'   1.5.
#' @param bp.lower Optional. Lower bound for the fitting of bp. Default is -10.
#' @param bp.upper Optional. Upper bound for the fitting of bp. Default is 15.
#' @param vBr.start Optional. Starting parameter for fitting of vBr. Default is
#'   0.05.
#' @param vBr.lower Optional. Lower bound for the fitting of vBr. Default is
#'   0.0001.
#' @param vBr.upper Optional. Upper bound for the fitting of vBr. Default is
#'   0.15.
#' @param vBt.start Optional. Starting parameter for fitting of vBt. Default is
#'   0.05.
#' @param vBt.lower Optional. Lower bound for the fitting of vBt. Default is
#'   0.0001.
#' @param vBt.upper Optional. Upper bound for the fitting of vBt. Default is
#'   0.15.
#' @param multstart_iter Number of iterations for starting parameters. Default
#'   is 1. For more information, see \code{\link[nls.multstart]{nls_multstart}}.
#'   If specified as 1 for any parameters, the original starting value will be
#'   used, and the multstart_lower and multstart_upper values ignored.
#' @param multstart_lower Optional. Lower bounds for starting parameters.
#'   Defaults to the lower bounds. Named list of whichever parameters' starting
#'   bounds should be altered.
#' @param multstart_upper Optional. Upper bounds for starting parameters.
#'   Defaults to the upper bounds. Named list of whichever parameters' starting
#'   bounds should be altered.
#' @param printvals Optional. This displays the parameter values for each
#'   iteration of the model. This is useful for debugging and changing starting
#'   values and upper and lower bounds for parameters.
#'
#' @return A list with a data frame of the fitted parameters \code{out$par},
#'   their percentage standard errors \code{out$par.se}, the model fit object
#'   \code{out$fit}, the model weights \code{out$weights}, and a dataframe
#'   containing the TACs both of the data and the fitted values \code{out$tacs}.
#'
#' @examples
#'
#' # Note: Reference region models, and irreversible binding models, should not
#' # be used for PBR28 - this is just to demonstrate function
#'
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' reftac <- pbr28$tacs[[2]]$CBL
#' roitac <- pbr28$tacs[[2]]$STR
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' input <- pbr28$input[[2]]
#'
#' newvals <- shift_timings(
#'   t_tac,
#'   roitac,
#'   input,
#'   inpshift = 0
#' )
#'
#' bloodtac <- pracma::interp1(newvals$input$Time, newvals$input$Blood, t_tac)
#'
#' fit <- srtm_v(t_tac, reftac, roitac, bloodtac, weights)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Tomasi, G., Edison, P., ... & Turkheimer, F. E. (2008). Novel
#'   reference region model reveals increased microglial and reduced vascular
#'   binding of 11C-(R)-PK11195 in patients with Alzheimer's disease. Journal of
#'   Nuclear Medicine, 49(8), 1249-1256.
#'
#' @export

srtm_v <- function(t_tac, reftac, roitac, bloodtac, weights = NULL, vBr = NULL, frameStartEnd = NULL,
                   R1.start = 1, R1.lower = 0.0001, R1.upper = 10,
                   k2.start = 0.1, k2.lower = 0.0001, k2.upper = 1,
                   bp.start = 1.5, bp.lower = -10, bp.upper = 15,
                   vBr.start = 0.05, vBr.lower = 0.0001, vBr.upper = 0.15,
                   vBt.start = 0.05, vBt.lower = 0.0001, vBt.upper = 0.15,
                   multstart_iter = 1, multstart_lower = NULL, multstart_upper = NULL,
                   printvals = F) {


  # Tidying

  tidyinput <- tidyinput_ref(t_tac, reftac, roitac, weights, frameStartEnd)
  tidyinput_blood <- tidyinput_ref(t_tac, reftac, bloodtac, weights, frameStartEnd)

  modeldata <- tidyinput
  modeldata$bloodtac <- tidyinput_blood$roitac


  # Parameters

  start <- c(R1 = R1.start, k2 = k2.start, bp = bp.start, vBr = vBr.start, vBt = vBt.start)
  lower <- c(R1 = R1.lower, k2 = k2.lower, bp = bp.lower, vBr = vBr.lower, vBt = vBt.lower)
  upper <- c(R1 = R1.upper, k2 = k2.upper, bp = bp.upper, vBr = vBr.upper, vBt = vBt.upper)

  vBr_fitted <- T
  if (!is.null(vBr)) {
    vBr_fitted <- F
    start[which(names(start) == "vBr")] <- vBr
    lower[which(names(lower) == "vBr")] <- vBr
    upper[which(names(upper) == "vBr")] <- vBr
  }

  multstart_pars <- fix_multstartpars(
    start, lower, upper, multstart_iter,
    multstart_lower, multstart_upper
  )
  multstart_upper <- multstart_pars$multstart_upper
  multstart_lower <- multstart_pars$multstart_lower



  # Solution

  if (prod(multstart_iter) == 1) {
    output <- minpack.lm::nlsLM(
      roitac ~ srtm_v_model(t_tac, reftac, bloodtac, R1, k2, bp, vBr, vBt),
      data = modeldata,
      start = start, lower = lower, upper = upper,
      weights = weights, control = minpack.lm::nls.lm.control(maxiter = 200),
      trace = printvals
    )
  } else {
    output <- nls.multstart::nls_multstart(
      roitac ~ srtm_v_model(t_tac, reftac, bloodtac, R1, k2, bp, vBr, vBt),
      data = modeldata, supp_errors = "Y",
      start_lower = multstart_lower,
      start_upper = multstart_upper,
      iter = multstart_iter, convergence_count = FALSE,
      lower = lower, upper = upper, modelweights = weights
    )
  }

  # Output

  tacs <- data.frame(
    Time = modeldata$t_tac, Reference = modeldata$reftac,
    Blood = modeldata$bloodtac, Target = modeldata$roitac,
    Target_fitted = as.numeric(fitted(output))
  )

  par <- as.data.frame(as.list(coef(output)))

  par.se <- par
  par.se[1,] <- purrr::map_dbl(names(coef(output)), ~ get_se(output, .x))
  names(par.se) <- paste0(names(par.se), ".se")

  out <- list(
    par = par, par.se = par.se,
    fit = output, weights = modeldata$weights, tacs = tacs,
    model = "srtm_v"
  )

  class(out) <- c("srtm_v", "kinfit")

  return(out)
}


#' Model: Simplified Reference Tissue Model with Blood Volumes
#'
#' This is the SRTM model itself by which predicted values are generated.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a zero.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame.
#' @param bloodtac Numeric vector of radioactivity concentrations in the blood for each frame.
#' @param R1 Parameter value for R1
#' @param k2 Parameter value for k2
#' @param bp Parameter value for bp
#' @param vBr Parameter value for vBr
#' @param vBt Parameter value for vBt
#'
#' @return A numeric vector of the predicted values of the TAC in the target region.
#'
#' @examples
#' \dontrun{
#' srtm_v_model(t_tac, reftac, R1 = 0.9, k2 = 0.1, bp = 1.5, vBr = 0.05, vBt = 0.05)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Tomasi, G., Edison, P., ... & Turkheimer, F. E. (2008). Novel
#'   reference region model reveals increased microglial and reduced vascular
#'   binding of 11C-(R)-PK11195 in patients with Alzheimer's disease. Journal of
#'   Nuclear Medicine, 49(8), 1249-1256.
#'
#' @export


srtm_v_model <- function(t_tac, reftac, bloodtac, R1, k2, bp, vBr, vBt) {
  interptime <- pracma::linspace(min(t_tac), max(t_tac), 1024)
  step <- interptime[2] - interptime[1]

  iref <- pracma::interp1(t_tac, reftac, interptime, method = "linear")
  iblood <- pracma::interp1(t_tac, bloodtac, interptime, method = "linear")

  a <- vBt * iblood
  b <- (1 - vBt) / (1 - vBr)

  c <- R1 * (iref - vBr * iblood)

  d <- k2 - (k2 * R1) / (1 + bp)

  e <- iref - vBr * iblood
  f <- exp((-k2 / (1 + bp)) * interptime)

  ef <- kinfit_convolve(e, f, step)

  itac <- a + b * (c + d * ef)

  outtac <- pracma::interp1(interptime, itac, t_tac)

  return(outtac)
}


#' Plot: Simplified Reference Tissue Model with Blood Volumes
#'
#' Function to visualise the fit of the SRTM_V model to data.
#'
#' @param srtmvout The output object of the SRTM fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#' @param refname Optional. The name of the Reference Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' # Note: Reference region models, and irreversible binding models, should not
#' # be used for PBR28 - this is just to demonstrate function
#'
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' reftac <- pbr28$tacs[[2]]$CBL
#' roitac <- pbr28$tacs[[2]]$STR
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' input <- pbr28$input[[2]]
#'
#' newvals <- shift_timings(
#'   t_tac,
#'   roitac,
#'   input,
#'   inpshift = 0
#' )
#'
#' bloodtac <- pracma::interp1(newvals$input$Time, newvals$input$Blood, t_tac)
#'
#' fit <- srtm_v(t_tac, reftac, roitac, bloodtac, weights)
#'
#' plot_srtm_vfit(fit)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_srtm_vfit <- function(srtmvout, roiname = NULL, refname = NULL) {
  measured <- data.frame(
    Time = srtmvout$tacs$Time,
    Reference = srtmvout$tacs$Reference,
    Blood = srtmvout$tacs$Blood,
    ROI.measured = srtmvout$tacs$Target,
    Weights = weights(srtmvout$fit)
  )

  fitted <- data.frame(
    Time = srtmvout$tacs$Time,
    ROI.fitted = srtmvout$tacs$Target_fitted,
    Weights = weights(srtmvout$fit)
  )

  if (is.null(roiname)) {
    roiname <- "ROI"
  }
  if (is.null(refname)) {
    refname <- "Reference"
  }

  measured <- plyr::rename(measured, c(
    "ROI.measured" = paste0(roiname, ".measured"),
    "Reference" = refname
  ))

  fitted <- plyr::rename(fitted, c("ROI.fitted" = paste0(roiname, ".fitted")))

  tidymeasured <- tidyr::gather(
    measured,
    key = Region, value = Radioactivity,
    -Time, -Weights, -Blood, factor_key = F
  )

  tidyfitted <- tidyr::gather(
    fitted,
    key = Region, value = Radioactivity,
    -Time, -Weights, factor_key = F
  )


  Region <- forcats::fct_inorder(factor(c(tidymeasured$Region, tidyfitted$Region)))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  outplot <- ggplot(tidymeasured, aes(x = Time, y = Radioactivity, colour = Region)) +
    geom_point(data = tidymeasured, aes(shape = "a", size = Weights)) +
    geom_line(data = tidyfitted) +
    geom_point(data = measured, aes(x = Time, y = Blood), colour = "black") +
    geom_line(data = measured, aes(x = Time, y = Blood), colour = "black") +
    guides(shape = FALSE, color = guide_legend(order = 1)) + colScale +
    scale_size(range = c(1, 3)) +
    coord_cartesian(ylim = c(0, max(tidymeasured$Radioactivity) * 1.5))

  return(outplot)
}
