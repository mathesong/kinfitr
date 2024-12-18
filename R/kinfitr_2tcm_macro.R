#' Two Tissue Compartment Model with Alternative Macroparameter Parameterisation
#'
#' Function to fit the Two Tissue Compartment Model with an alternative
#' parameterisation to data. This gives some more fine-grained control over
#' setting parameter limits and setting parameter values from previously
#' estimated quantities.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the
#'   time halfway through the frame as well as a zero. If a time zero frame is
#'   not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target
#'   tissue for each frame. We include zero at time zero: if not included, it is
#'   added.
#' @param input Data frame containing the blood, plasma, and parent fraction
#'   concentrations over time.  This can be generated using the
#'   \code{blood_interp} function.
#' @param weights Optional. Numeric vector of the weights assigned to each frame
#'   in the fitting. We include zero at time zero: if not included, it is added.
#'   If not specified, uniform weights will be used.
#' @param inpshift Optional. The number of minutes by which to shift the timing
#'   of the input data frame forwards or backwards. If not specified, this will
#'   be fitted, however this takes longer to compute. Recommended to perform
#'   once on a large ROI for each measurement, and to specify this value for the
#'   remainder of the regions. Or perhaps even to use a simpler model, such as
#'   1TCM.
#' @param vB  Optional. The blood volume fraction.  If not specified, this will
#'   be fitted. Recommended to perform once on a large ROI for each measurement,
#'   and to specify this value for the remainder of the regions.
#' @param Vnd Optional. This parameter can be used to fix the non-displaceable
#'   volume of distribution, i.e. the K1/k2 ratio to a specific value.  If not
#'   specified, this will be fitted.
#' @param BPp Optional. This parameter can be used to fix the BPp value, i.e.
#'   K1k3/k2k4 to a specific value.  If not specified, this will be fitted.
#' @param k4 Optional. This parameter can be used to fix k4 to a specific value.
#'   If not specified, this will be fitted.
#' @param frameStartEnd Optional: This allows one to specify the beginning and
#'   final frame to use for modelling, e.g. c(1,20). This is to assess time
#'   stability.
#' @param K1.start Optional. Starting parameter for fitting of K1. Default is
#'   0.1.
#' @param K1.lower Optional. Lower bound for the fitting of K1. Default is
#'   0.0001.
#' @param K1.upper Optional. Upper bound for the fitting of K1. Default is 1.
#' @param Vnd.start Optional. Starting parameter for fitting of Vnd. Default is
#'   1.
#' @param Vnd.lower Optional. Lower bound for the fitting of Vnd. Default is
#'   0.0001.
#' @param Vnd.upper Optional. Upper bound for the fitting of Vnd. Default is 10.
#' @param BPp.start Optional. Starting parameter for fitting of BPp. Default is
#'   1.
#' @param BPp.lower Optional. Lower bound for the fitting of BPp. Default is
#'   0.0001.
#' @param BPp.upper Optional. Upper bound for the fitting of BPp. Default is 50.
#' @param k4.start Optional. Starting parameter for fitting of k4. Default is
#'   0.1.
#' @param k4.lower Optional. Lower bound for the fitting of k4. Default is
#'   0.0001.
#' @param k4.upper Optional. Upper bound for the fitting of k4. Default is 0.5.
#' @param inpshift.start Optional. Starting parameter for fitting of inpshift.
#'   Default is 0.
#' @param inpshift.lower Optional. Lower bound for the fitting of inpshift.
#'   Default is -0.5.
#' @param inpshift.upper Optional. Upper bound for the fitting of inpshift.
#'   Default is 0.5.
#' @param vB.start Optional. Starting parameter for fitting of vB. Default is
#'   0.05.
#' @param vB.lower Optional. Lower bound for the fitting of vB. Default is 0.01.
#' @param vB.upper Optional. Upper bound for the fitting of vB. Default is 0.1.
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
#'   their percentage standard errors (scaled so that 1 represents 100\%)
#'   \code{out$par.se}, the model fit object \code{out$fit}, a dataframe
#'   containing the TACs both of the data and the fitted values \code{out$tacs},
#'   the blood input data frame after time shifting \code{input}, a vector of
#'   the weights \code{out$weights}, a logical of whether the inpshift was
#'   fitted \code{inpshift_fitted} and a logical of whether the vB was fitted
#'   \code{vB}.
#'
#' @examples
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' tac <- pbr28$tacs[[2]]$FC
#' weights <- pbr28$tacs[[2]]$Weights
#'
#' input <- blood_interp(
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cbl_dispcorr,
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cpl_metabcorr,
#'   t_parentfrac = 1, parentfrac = 1
#' )
#'
#' fit1 <- twotcm(t_tac, tac, input, weights)
#' fit2 <- twotcm(t_tac, tac, input, weights, inpshift = 0.1, vB = 0.05)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

twotcm_macro <- function(t_tac, tac, input, weights = NULL, inpshift = NULL,
                         vB = NULL, Vnd = NULL, BPp = NULL, k4 = NULL,
                         frameStartEnd = NULL,
                         K1.start = 0.1, K1.lower = 0.0001, K1.upper = 1,
                         Vnd.start = 1, Vnd.lower = 0.0001, Vnd.upper = 10,
                         BPp.start = 1, BPp.lower = 0.0001, BPp.upper = 50,
                         k4.start = 0.1, k4.lower = 0.0001, k4.upper = 0.5,
                         inpshift.start = 0, inpshift.lower = -0.5, inpshift.upper = 0.5,
                         vB.start = 0.05, vB.lower = 0.01, vB.upper = 0.1,
                         multstart_iter = 1, multstart_lower = NULL, multstart_upper = NULL,
                         printvals = F) {

  # Tidying

  tidyinput <- tidyinput_art(t_tac, tac, weights, frameStartEnd)

  modeldata <- list(
    t_tac = tidyinput$t_tac,
    tac = tidyinput$tac,
    weights = tidyinput$weights,
    input = input
  )

  # Parameters

  start <- c(
    K1 = K1.start, Vnd = Vnd.start, BPp = BPp.start, k4 = k4.start,
    inpshift = inpshift.start, vB = vB.start
  )
  lower <- c(
    K1 = K1.lower, Vnd = Vnd.lower, BPp = BPp.lower, k4 = k4.lower,
    inpshift = inpshift.lower, vB = vB.lower
  )
  upper <- c(
    K1 = K1.upper, Vnd = Vnd.upper, BPp = BPp.upper, k4 = k4.upper,
    inpshift = inpshift.upper, vB = vB.upper
  )

  multstart_pars <- fix_multstartpars(
    start, lower, upper, multstart_iter,
    multstart_lower, multstart_upper
  )
  multstart_upper <- multstart_pars$multstart_upper
  multstart_lower <- multstart_pars$multstart_lower

  # Sort out whether vB is to be fitted
  vB_fitted <- T
  if (!is.null(vB)) {
    vB_fitted <- F

    par_keepindex <- names(start) != "vB"
    start <- start[par_keepindex]
    lower <- lower[par_keepindex]
    upper <- upper[par_keepindex]
    multstart_lower <- multstart_lower[par_keepindex]
    multstart_upper <- multstart_upper[par_keepindex]
  }

  # Sort out whether Vnd is to be fitted
  Vnd_fitted <- T
  if (!is.null(Vnd)) {
    Vnd_fitted <- F

    par_keepindex <- names(start) != "Vnd"
    start <- start[par_keepindex]
    lower <- lower[par_keepindex]
    upper <- upper[par_keepindex]
    multstart_lower <- multstart_lower[par_keepindex]
    multstart_upper <- multstart_upper[par_keepindex]
  }

  # Sort out whether BPp is to be fitted
  BPp_fitted <- T
  if (!is.null(BPp)) {
    BPp_fitted <- F

    par_keepindex <- names(start) != "BPp"
    start <- start[par_keepindex]
    lower <- lower[par_keepindex]
    upper <- upper[par_keepindex]
    multstart_lower <- multstart_lower[par_keepindex]
    multstart_upper <- multstart_upper[par_keepindex]
  }

  # Sort out whether k4 is to be fitted
  k4_fitted <- T
  if (!is.null(k4)) {
    k4_fitted <- F

    par_keepindex <- names(start) != "k4"
    start <- start[par_keepindex]
    lower <- lower[par_keepindex]
    upper <- upper[par_keepindex]
    multstart_lower <- multstart_lower[par_keepindex]
    multstart_upper <- multstart_upper[par_keepindex]
  }


  # Solution - Delay Already Fitted

  if (!is.null(inpshift)) {
    inpshift_fitted <- F

    par_keepindex <- names(start) != "inpshift"
    start <- start[par_keepindex]
    lower <- lower[par_keepindex]
    upper <- upper[par_keepindex]
    multstart_lower <- multstart_lower[par_keepindex]
    multstart_upper <- multstart_upper[par_keepindex]

    newvals <- shift_timings(
      modeldata$t_tac,
      modeldata$tac,
      modeldata$input,
      inpshift
    )

    modeldata$t_tac <- newvals$t_tac
    modeldata$tac <- newvals$tac
    modeldata$input <- newvals$input

    formula <- paste(
      "tac ~ twotcm_macro_model(t_tac, input, K1",
      ", Vnd", ifelse( Vnd_fitted, yes = "", no = paste0("=", Vnd)),
      ", BPp",ifelse( BPp_fitted, yes = "", no = paste0("=", BPp)),
      ", k4",ifelse( k4_fitted, yes = "", no = paste0("=", k4)),
      ", vB", ifelse( vB_fitted, yes = "", no = paste0("=", vB)),
      ")", sep=""
    )

    if (prod(multstart_iter) == 1) {
      output <- minpack.lm::nlsLM(
        as.formula(formula),
        data = modeldata,
        start = start, lower = lower, upper = upper,
        weights = weights, control = minpack.lm::nls.lm.control(maxiter = 200),
        trace = printvals
      )
    } else {
      output <- nls.multstart::nls_multstart(
        as.formula(formula),
        data = modeldata,
        supp_errors = "Y",
        start_lower = multstart_lower,
        start_upper = multstart_upper,
        iter = multstart_iter, convergence_count = FALSE,
        lower = lower, upper = upper, modelweights = weights
      )
    }
  }

  # Solution - Fitting the Delay

  if (is.null(inpshift)) {
    inpshift_fitted <- T

    formula <- paste(
      "tac ~ twotcm_macro_fitDelay_model(t_tac, input, K1",
      ", Vnd", ifelse( Vnd_fitted, yes = "", no = paste0("=", Vnd)),
      ", BPp",ifelse( BPp_fitted, yes = "", no = paste0("=", BPp)),
      ", k4",ifelse( k4_fitted, yes = "", no = paste0("=", k4)),
      ", inpshift",
      ", vB", ifelse( vB_fitted, yes = "", no = paste0("=", vB)),
      ")", sep=""
    )

    if (prod(multstart_iter) == 1) {
      output <- minpack.lm::nlsLM(
        as.formula(formula),
        data = modeldata,
        start = start, lower = lower, upper = upper,
        weights = weights, control = minpack.lm::nls.lm.control(maxiter = 200),
        trace = printvals
      )
    } else {
      output <- nls.multstart::nls_multstart(
        as.formula(formula),
        data = modeldata,
        supp_errors = "Y",
        start_lower = multstart_lower,
        start_upper = multstart_upper,
        iter = multstart_iter, convergence_count = FALSE,
        lower = lower, upper = upper, modelweights = weights
      )
    }
  }

  # Check for parameters hitting limits

  limcheck_u <- purrr::map2_lgl(round(upper,3), round(coef(output),3), identical)
  limcheck_l <- purrr::map2_lgl(round(lower,3), round(coef(output),3), identical)
  limcheck <- limcheck_u + limcheck_l
  limcheck <- limcheck==1

  if(
    any(limcheck)
  ) {
    warning(
      paste0(
        "\nFitted parameters are hitting upper or lower limit bounds. Consider \n",
        "either modifying the upper and lower limit boundaries, or else using \n",
        "multstart when fitting the model (see the function documentation).\n") )
  }

  # Output

  if (inpshift_fitted == T) {
    newvals <- shift_timings(
      modeldata$t_tac,
      modeldata$tac, modeldata$input,
      as.numeric(coef(output)[["inpshift"]])
    )

    tacs <- data.frame(Time = newvals$t_tac, Target = newvals$tac, Target_fitted = as.numeric(fitted(output)))
    input <- newvals$input
  } else {
    tacs <- data.frame(Time = newvals$t_tac, Target = newvals$tac, Target_fitted = as.numeric(fitted(output)))
    input <- newvals$input
  }

  par <- as.data.frame(as.list(coef(output)))

  if (inpshift_fitted == F) par$inpshift <- inpshift
  if (Vnd_fitted == F)       par$Vnd <- Vnd
  if (BPp_fitted == F)       par$BPp <- BPp
  if (k4_fitted == F)        par$k4 <- k4
  if (vB_fitted == F)        par$vB <- vB

  par.se <- par
  par.se[1,] <- purrr::map_dbl(names(par), ~ get_se(output, .x))
  names(par.se) <- paste0(names(par.se), ".se")

  par$Vt <- (par$Vnd + par$BPp)
  par.se$Vt.se <- get_se(output, "Vnd + BPp")

  par$BPnd <- (par$BPp / par$Vnd)
  par.se$BPnd.se <- get_se(output, "BPp / Vnd")

  par$k3 <- (par$BPp / par$Vnd) * par$k4
  par.se$k3.se <- get_se(output, "(BPp / Vnd) * k4")

  par$k2 <- par$K1/par$Vnd
  par.se$k2.se <- get_se(output, "K1 / Vnd")



  out <- list(
    par = par, par.se = par.se,
    fit = output, tacs = tacs, input = input, weights = modeldata$weights,
    inpshift_fitted = inpshift_fitted,
    Vnd_fitted = Vnd_fitted, BPp_fitted = BPp_fitted,
    k4_fitted = k4_fitted, vB_fitted = vB_fitted,
    model = "2tcm"
  )

  class(out) <- c("2tcm", "kinfit")

  return(out)
}


#' Model: Two Tissue Compartment Model with a Macroparameter Parameterisation
#'
#' This is the Two Tissue Compartment Model model itself by which predicted values are generated.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a zero.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param K1 Parameter value for K1
#' @param Vnd Parameter value for Vnd
#' @param BPp Parameter value for BPp
#' @param k4 Parameter value for k4
#' @param vB Parameter value for vB
#'
#' @return A numeric vector of the predicted values of the TAC in the target region.
#'
#' @examples
#' \dontrun{
#' twotcm_macro_model(t_tac, input, K1 = 0.1, Vnd = 1, BPp = 2, k4 = 0.02, vB = 0.05)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

twotcm_macro_model <- function(t_tac, input, K1, Vnd, BPp, k4, vB) {
  interptime <- input$Time
  step <- interptime[2] - interptime[1]

  i_blood <- input$Blood
  aif <- input$AIF

  k2 <- K1 / Vnd
  k3 <- (BPp / Vnd) * k4

  delta <- sqrt((k2 + k3 + k4)^2 - 4 * k2 * k4)
  th1 <- (k2 + k3 + k4 + delta) / 2
  th2 <- (k2 + k3 + k4 - delta) / 2
  ph1 <- K1 * (th1 - k3 - k4) / delta
  ph2 <- K1 * (th2 - k3 - k4) / (-delta)

  a <- ph1 * exp(-th1 * interptime) + ph2 * exp(-th2 * interptime)
  b <- aif

  i_outtac <- kinfit_convolve(a, b, step)

  # Correction for vB
  i_outtac <- i_outtac * (1 - vB) + vB * i_blood

  outtac <- pracma::interp1(interptime, i_outtac, t_tac)

  return(outtac)
}


#' Model: Two Tissue Compartment Model with Delay with a Macroparameter Parameterisation
#'
#' This is the Two Tissue Compartment Model model itself by which predicted values are generated, which includes fitting of the
#' delay, inpshift.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a zero.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param K1 Parameter value for K1
#' @param Vnd Parameter value for Vnd
#' @param BPp Parameter value for BPp
#' @param k4 Parameter value for k4
#' @param inpshift Parameter value for inpshift, the delay.
#' @param vB Parameter value for vB
#'
#' @return A numeric vector of the predicted values of the TAC in the target region.
#'
#' @examples
#' \dontrun{
#' twotcm_macro_fitDelay_model(t_tac, input, K1 = 0.1, Vnd = 1, BPp = 2,
#'                             k4 = 0.02, inpshift = 0, vB = 0.05)
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

twotcm_macro_fitDelay_model <- function(t_tac, input, K1, Vnd, BPp, k4, inpshift, vB) {
  newvals <- shift_timings(t_tac, rep(1, length(t_tac)), input, inpshift) # Using ones instead of tac as don't need it

  t_tac <- newvals$t_tac

  i_blood <- newvals$input$Blood
  aif <- newvals$input$AIF

  interptime <- newvals$input$Time
  step <- interptime[2] - interptime[1]

  k2 <- K1 / Vnd
  k3 <- (BPp / Vnd) * k4

  delta <- sqrt((k2 + k3 + k4)^2 - 4 * k2 * k4)
  th1 <- (k2 + k3 + k4 + delta) / 2
  th2 <- (k2 + k3 + k4 - delta) / 2
  ph1 <- K1 * (th1 - k3 - k4) / delta
  ph2 <- K1 * (th2 - k3 - k4) / (-delta)

  a <- ph1 * exp(-th1 * interptime) + ph2 * exp(-th2 * interptime)
  b <- aif

  i_outtac <- kinfit_convolve(a, b, step)

  # Correction for vB
  i_outtac <- i_outtac * (1 - vB) + vB * i_blood

  outtac <- pracma::interp1(interptime, i_outtac, t_tac)

  return(outtac)
}
