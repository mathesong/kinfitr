#' Two Tissue Compartment Model with Irreversible Trapping
#'
#' Function to fit the Two Tissue Compartment Model to data. An irreversible model can also be specified by setting the
#' starting value, upper and lower bounds of k4 to 0.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero:
#' if not included, it is added. If not specified, uniform weights will be used.
#' @param inpshift Optional. The number of minutes by which to shift the timing of the input data frame forwards or backwards.
#' If not specified, this will be fitted, however this takes longer to compute. Recommended to perform once on a large ROI for
#' each measurement, and to specify this value for the remainder of the regions. Or perhaps even to use a simpler model, such as
#' 1TCM.
#' @param vB_fixed Optional. The blood volume fraction.  If not specified, this will be fitted. Recommended to perform once on a large ROI for
#' each measurement, and to specify this value for the remainder of the regions.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20).
#' This is to assess time stability.
#' @param K1.start Optional. Starting parameter for fitting of K1. Default is 0.1.
#' @param K1.lower Optional. Lower bound for the fitting of K1. Default is 0.0001.
#' @param K1.upper Optional. Upper bound for the fitting of K1. Default is 0.5.
#' @param k2.start Optional. Starting parameter for fitting of k2. Default is 0.1.
#' @param k2.lower Optional. Lower bound for the fitting of k2. Default is 0.0001.
#' @param k2.upper Optional. Upper bound for the fitting of k2. Default is 0.5.
#' @param k3.start Optional. Starting parameter for fitting of k3. Default is 0.1.
#' @param k3.lower Optional. Lower bound for the fitting of k3. Default is 0.0001.
#' @param k3.upper Optional. Upper bound for the fitting of k3. Default is 0.5.
#' @param k4.start Optional. Starting parameter for fitting of k4. Default is 0.1.
#' @param k4.lower Optional. Lower bound for the fitting of k4. Default is 0.0001.
#' @param k4.upper Optional. Upper bound for the fitting of k4. Default is 0.5.
#' @param Kb.start Optional. Starting parameter for fitting of Kb. Default is 0.25.
#' @param Kb.lower Optional. Lower bound for the fitting of Kb. Default is 0.0001.
#' @param Kb.upper Optional. Upper bound for the fitting of Kb. Default is 1.
#' @param inpshift.start Optional. Starting parameter for fitting of inpshift. Default is 0.
#' @param inpshift.lower Optional. Lower bound for the fitting of inpshift. Default is -0.5.
#' @param inpshift.upper Optional. Upper bound for the fitting of inpshift. Default is 0.5.
#' @param vB.start Optional. Starting parameter for fitting of vB. Default is 0.05.
#' @param vB.lower Optional. Lower bound for the fitting of vB. Default is 0.01.
#' @param vB.upper Optional. Upper bound for the fitting of vB. Default is 0.1.
#' @param multstart_iter Number of iterations for starting parameters. Default is 1.
#'   For more information, see \code{\link[nls.multstart]{nls_multstart}}. If
#'   specified as 1 for any parameters, the original starting value will be
#'   used, and the multstart_lower and multstart_upper values ignored.
#' @param multstart_lower Optional. Lower bounds for starting parameters. Defaults
#'   to the lower bounds.
#' @param multstart_upper Optional. Upper bounds for starting parameters. Defaults
#'   to the upper bounds.
#' @param printvals Optional. This displays the parameter values for each iteration of the
#' model. This is useful for debugging and changing starting values and upper and lower
#' bounds for parameters.
#'
#' @return A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit},
#' a dataframe containing the TACs both of the data and the fitted values \code{out$tacs},
#' the blood input data frame after time shifting \code{input}, a vector of the weights \code{out$weights},
#' a logical of whether the inpshift was fitted \code{inpshift_fitted} and a logical of whether the vB was
#' fitted \code{vB_fixed}.
#'
#' @examples
#' twotcm1k(t_tac, tac, input, weights=weights)
#' twotcm1k(t_tac, tac, input, weights=weights, inpshift=0.1, vB_fixed=0.05)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Rizzo, G., Veronese, M., Tonietto, M., Zanotti-Fregonara, P., Turkheimer, F. E., & Bertoldo, A. (2014). Kinetic modeling without accounting for the vascular component impairs the quantification of [11C] PBR28 brain PET data. Journal of Cerebral Blood Flow & Metabolism, 34(6), 1060-1069.
#'
#' @export

twotcm1k <- function(t_tac, tac, input, weights, inpshift, vB_fixed, frameStartEnd,
                   K1.start = 0.1 , K1.lower = 0.0001 , K1.upper = 0.5 ,
                   k2.start = 0.1 , k2.lower = 0.0001 , k2.upper = 0.5 ,
                   k3.start = 0.1 , k3.lower = 0.0001 , k3.upper = 0.5 ,
                   k4.start = 0.1 , k4.lower = 0.0001 , k4.upper = 0.5 ,
                   Kb.start = 0.25 , Kb.lower = 0.0001 , Kb.upper = 1 ,
                   inpshift.start = 0 , inpshift.lower= -0.5 , inpshift.upper = 0.5 ,
                   vB.start = 0.05 , vB.lower = 0.01 , vB.upper = 0.1,
                   multstart_iter=1, multstart_lower, multstart_upper,
                   printvals=F) {

  # Tidying

  tidyinput <- tidyinput_art(t_tac, tac, weights, frameStartEnd)

  modeldata <- list(t_tac = tidyinput$t_tac,
                    tac = tidyinput$tac,
                    weights = tidyinput$weights,
                    input = input)

  # Parameters

  start <- c(K1 = K1.start, k2 = k2.start, k3 = k3.start, k4 = k4.start,
             Kb = Kb.start, inpshift = inpshift.start, vB = vB.start)
  lower <- c(K1 = K1.lower, k2 = k2.lower, k3 = k3.lower, k4 = k4.lower,
             Kb = Kb.lower, inpshift = inpshift.lower, vB = vB.lower)
  upper <- c(K1 = K1.upper, k2 = k2.upper, k3 = k3.upper, k4 = k4.upper,
             Kb = Kb.upper, inpshift = inpshift.upper, vB = vB.upper)

  vB_fitted = T
  if(!missing(vB_fixed)) {
    vB_fitted = F

    start[which(names(start)=='vB')] <- vB_fixed
    lower[which(names(lower)=='vB')] <- vB_fixed
    upper[which(names(upper)=='vB')] <- vB_fixed
  }

  if( length(multstart_iter) == length(start) || length(multstart_iter) == 1 ) {

    ### Missing multstart boundaries
    if(missing(multstart_lower)) {
      multstart_lower = lower
    }
    if(missing(multstart_upper)) {
      multstart_upper = upper
    }

    ### No multstart for some variables ###
    if( any(multstart_iter==1) ) {
      non_iterable <- which(multstart_iter==1)
      multstart_lower[non_iterable] = start[non_iterable]
      multstart_upper[non_iterable] = start[non_iterable]
    }

    ### Incorrect multstart boundaries
    if( length(multstart_lower) != length(start) || length(multstart_upper) != length(start) ) {
      stop('multstart_lower and multstart_upper should be the same length as the number of parameters')
    }
  } else {
    stop('multstart_iter should be of length 1 or of the same length as the number of parameters')
  }


  # Solution - Delay Already Fitted

  if(!missing(inpshift)) {

    inpshift_fitted = F
    multstart_lower <- multstart_lower[names(multstart_lower)!='inpshift']
    multstart_upper <- multstart_upper[names(multstart_upper)!='inpshift']

    newvals <- shift_timings(t_tac, tac, input, inpshift)

    modeldata$t_tac <- newvals$t_tac
    modeldata$tac <- newvals$tac
    modeldata$input <- newvals$input

    start <- start[-which(names(start)=='inpshift')]
    lower <- lower[-which(names(lower)=='inpshift')]
    upper <- upper[-which(names(upper)=='inpshift')]

    if( prod(multstart_iter) == 1 ) {

      output <- minpack.lm::nlsLM(tac ~ twotcm1k_model(t_tac, input, K1, k2, k3, k4, Kb, vB),
                                  data=modeldata,
                                  start =  start, lower = lower, upper = upper,
                                  weights = weights, control = minpack.lm::nls.lm.control(maxiter = 200),
                                  trace = printvals)

    } else {

      output <- nls.multstart::nls_multstart(tac ~ twotcm1k_model(t_tac, input, K1, k2, k3, k4, Kb, vB),
                                             data=modeldata,
                                             supp_errors = 'Y',
                                             start_lower = multstart_lower,
                                             start_upper = multstart_upper,
                                             iter = multstart_iter, convergence_count = FALSE,
                                             modelweights=weights)

    }
  }

  # Solution - Fitting the Delay

  if(missing(inpshift)) {

    inpshift_fitted = T

    if( prod(multstart_iter) == 1 ) {

      output <- minpack.lm::nlsLM(tac ~ twotcm1k_fitDelay_model(t_tac, input, K1, k2, k3, k4, Kb, inpshift, vB),
                                  data=modeldata,
                                  start =  start, lower = lower, upper = upper,
                                  weights = weights, control = minpack.lm::nls.lm.control(maxiter = 200),
                                  trace = printvals)

    } else {

      output <- nls.multstart::nls_multstart(tac ~ twotcm1k_fitDelay_model(t_tac, input, K1, k2, k3, k4, Kb, inpshift, vB),
                                             data=modeldata,
                                             supp_errors = 'Y',
                                             start_lower = multstart_lower,
                                             start_upper = multstart_upper,
                                             iter = multstart_iter, convergence_count = FALSE,
                                             modelweights=weights)

    }
  }


  # Output

  if(inpshift_fitted == T) {
    newvals <- shift_timings(t_tac, tac, input, as.numeric(coef(output)[['inpshift']]))
    tacs <- data.frame(Time =newvals$t_tac, Target = newvals$tac, Target_fitted=as.numeric(fitted(output)))
    input <- newvals$input
  } else {
    tacs <- data.frame(Time =newvals$t_tac, Target = newvals$tac, Target_fitted=as.numeric(fitted(output)))
    input <- newvals$input
  }

  par = as.data.frame(as.list(coef(output)))

  if(inpshift_fitted == F) par$inpshift=inpshift

  par.se = as.data.frame(as.list(sqrt(abs(vcov(output)[,1]))))
  names(par.se) = paste0(names(par.se), '.se')

  par$Vt = (par$K1 / par$k2) * (1 + par$k3/par$k4)

  par.se$Vt.se = par$Vt * (par.se$K1.se/par$K1 + par.se$k2.se/par$k2) +
    par$Vt * (par.se$K1.se/par$K1 + par.se$k2.se/par$k2 + par.se$k3.se/par$k3 + par.se$k4.se/par$k4)



  out <- list(par = par, par.se = par.se,
              fit = output, tacs = tacs, input = input, weights = weights,
              inpshift_fitted = inpshift_fitted, vB_fitted = vB_fitted, model='2tcm1k')

  return(out)
}


#' Model: Two Tissue Compartment Model with Irreversible Trapping
#'
#' This is the Two Tissue Compartment Model with Irreversible Trapping model itself by which predicted values are generated.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a zero.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param K1 Parameter value for K1
#' @param k2 Parameter value for k2
#' @param k3 Parameter value for k3
#' @param k4 Parameter value for k4
#' @param Kb Parameter value for Kb
#' @param vB Parameter value for vB
#'
#' @return A numeric vector of the predicted values of the TAC in the target region.
#'
#' @examples
#' twotcm1k_model(t_tac, input, K1=0.1, k2=0.08, k3=0.05, k4=0.02, Kb=0.25, vB=0.05)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Rizzo, G., Veronese, M., Tonietto, M., Zanotti-Fregonara, P., Turkheimer, F. E., & Bertoldo, A. (2014). Kinetic modeling without accounting for the vascular component impairs the quantification of [11C] PBR28 brain PET data. Journal of Cerebral Blood Flow & Metabolism, 34(6), 1060-1069.
#'
#' @export

twotcm1k_model <- function(t_tac, input, K1, k2, k3, k4, Kb, vB) {

  interptime <- input$Time
  step <- interptime[2] - interptime[1]

  t_inp <- interptime
  i_blood <- input$blood
  i_plasma <- input$plasma
  i_parentfrac <- input$parentfrac

  i_inp <- i_plasma * i_parentfrac

  alpha <- ( (k2 + k3+ k4) - sqrt( (k2+k3+k4)^2 - (4*k2*k4)) ) / 2

  beta <- ( (k2 + k3+ k4) + sqrt( (k2+k3+k4)^2 - (4*k2*k4)) ) / 2


  A <- ( (k3+k4-alpha) / (beta-alpha) ) * exp(-alpha*interptime)

  B <- ( (beta - k3 - k4) / (beta-alpha) ) * exp(-beta*interptime)

  C <- vB * Kb * as.numeric(pracma::cumtrapz(interptime, i_inp))

  D <- vB * i_blood


  A_conv <- kinfit_convolve(A, i_inp, step)

  B_conv <- kinfit_convolve(B, i_inp, step)



  i_outtac <- (1-vB) * ( K1 * (A_conv + B_conv) ) + C + D

  outtac <- pracma::interp1(interptime, i_outtac, t_tac)

  return(outtac)
}


#' Model: Two Tissue Compartment Model with Irreversible Trapping with Delay Fitting
#'
#' This is the Two Tissue Compartment Model with Irreversible Trapping model itself by which predicted values are generated, which includes fitting of the
#' delay, inpshift.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a zero.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param K1 Parameter value for K1
#' @param k2 Parameter value for k2
#' @param k3 Parameter value for k3
#' @param k4 Parameter value for k4
#' @param Kb Parameter value for Kb
#' @param inpshift Parameter value for inpshift, the delay.
#' @param vB Parameter value for vB
#'
#' @return A numeric vector of the predicted values of the TAC in the target region.
#'
#' @examples
#' twotcm1k_fitDelay_model(t_tac, input, K1=0.1, k2=0.08, k3=0.05, k4=0.02, Kb=0.25, inpshift=0.1, vB=0.05)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Rizzo, G., Veronese, M., Tonietto, M., Zanotti-Fregonara, P., Turkheimer, F. E., & Bertoldo, A. (2014). Kinetic modeling without accounting for the vascular component impairs the quantification of [11C] PBR28 brain PET data. Journal of Cerebral Blood Flow & Metabolism, 34(6), 1060-1069.
#'
#' @export

twotcm1k_fitDelay_model <- function(t_tac, input, K1, k2, k3, k4, Kb, inpshift, vB) {

  newvals <- shift_timings(t_tac, rep(1,length(t_tac)), input, inpshift)

  t_tac <- newvals$t_tac

  t_inp <- newvals$input$Time
  i_blood <- newvals$input$blood
  i_plasma <- newvals$input$plasma
  i_parentfrac <- newvals$input$parentfrac

  interptime <- newvals$input$Time
  step <- interptime[2] - interptime[1]

  i_inp <- i_plasma * i_parentfrac

  alpha <- ( (k2 + k3+ k4) - sqrt( (k2+k3+k4)^2 - (4*k2*k4)) ) / 2

  beta <- ( (k2 + k3+ k4) + sqrt( (k2+k3+k4)^2 - (4*k2*k4)) ) / 2


  A <- ( (k3+k4-alpha) / (beta-alpha) ) * exp(-alpha*interptime)

  B <- ( (beta - k3 - k4) / (beta-alpha) ) * exp(-beta*interptime)

  C <- vB * Kb * as.numeric(pracma::cumtrapz(interptime, i_inp))

  D <- vB * i_blood


  A_conv <- kinfit_convolve(A, i_inp, step)

  B_conv <- kinfit_convolve(B, i_inp, step)



  i_outtac <- (1-vB) * ( K1 * (A_conv + B_conv) ) + C + D

  outtac <- pracma::interp1(interptime, i_outtac, t_tac)

  return(outtac)
}


#' Plot: Two Tissue Compartment Model with Irreversible Trapping
#'
#' Function to visualise the fit of the Two Tissue Compartment Model with Irreversible Trapping to data.
#'
#' @param twotcm1kout The output object of the 2TCM1k fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' plot_2tcm1kfit(twotcm1kout)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @references Rizzo, G., Veronese, M., Tonietto, M., Zanotti-Fregonara, P., Turkheimer, F. E., & Bertoldo, A. (2014). Kinetic modeling without accounting for the vascular component impairs the quantification of [11C] PBR28 brain PET data. Journal of Cerebral Blood Flow & Metabolism, 34(6), 1060-1069.
#'
#' @import ggplot2
#'
#' @export

plot_2tcm1kfit <- function(twotcm1kout, roiname) {

  if(missing(roiname)) {roiname = 'ROI'}

  measureddf <- data.frame(Time = twotcm1kout$tacs$Time,
                           Radioactivity = twotcm1kout$tacs$Target,
                           Weights = weights(twotcm1kout$fit),
                           Region=paste0(roiname, '.Measured'))

  inputdf <- data.frame(Time = twotcm1kout$input$Time,
                        Radioactivity = twotcm1kout$input$plasma*twotcm1kout$input$parentfrac,
                        Weights = 1,
                        Region='AIF')

  i_fit <- predict(twotcm1kout$fit, newdata = list(t_tac = twotcm1kout$input$Time,
                                                   tac = pracma::interp1(twotcm1kout$tacs$Time,
                                                                         twotcm1kout$tacs$Target,
                                                                         twotcm1kout$input$Time,
                                                                         method = "linear")))


  fitdf <- data.frame(Time = twotcm1kout$input$Time, Radioactivity = i_fit,
                      Weights = 1, Region = paste0(roiname, ".Fitted"))

  plotdf <- rbind(inputdf, measureddf, fitdf)
  plotdf$Region <- forcats::fct_inorder(factor(plotdf$Region))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(plotdf$Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  outplot = ggplot(plotdf, aes(x=Time, y=Radioactivity, colour=Region)) + colScale +
    geom_point(data=subset(plotdf, plotdf$Region == paste0(roiname, '.Measured')), aes(shape='a', size=Weights)) +
    geom_line(data=subset(plotdf, plotdf$Region != paste0(roiname, '.Measured'))) +
    guides(shape=FALSE, color=guide_legend(order=1)) + scale_size(range=c(1,3)) + coord_cartesian(ylim=c(0,max(measureddf$Radioactivity)*1.5))

  #print(outplot)
  return(outplot)
}
