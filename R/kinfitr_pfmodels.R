#' Hill Function Model for Parent Fraction
#'
#' This is the model function for fitting of the classic Hill function.
#'
#' @param time Time in seconds.
#' @param a Parameter A.
#' @param b Parameter B.
#' @param c Parameter C as the log10 of the 'classic' c parameter.
#' @param ppf0 The starting point of the parent fraction curve.
#' @param delay The delay of the metabolism curve.
#'
#' @return Predicted values
#' @export
#'
#' @examples
#' metab_hill_model(seq(0, 60 * 60, by = 120), 0.05, 2.6, 6.9, 1, 0)
metab_hill_model <- function(time, a, b, c, ppf0 = 1, delay = 0) {
  tcorr <- time - delay
  t_before <- tcorr[ which(!(tcorr > 0)) ]
  t_after <- tcorr[ which(tcorr > 0) ]

  ind1_out <- rep(ppf0, length.out = length(t_before))
  ind2_out <- ppf0 - (
    ((ppf0 - a) * t_after^b) /
      (10^c + (t_after)^b)
  )

  out <- c(ind1_out, ind2_out)

  return(out)
}



#' Fit a Hill Function to Model Parent Fraction
#'
#' This function fits a classic Hill Function to parent fraction data.
#'
#' @param time Time in seconds.
#' @param parentFraction Measured values of parent fraction.
#' @param fit_ppf0 Should the starting plasma parent fraction be fitted? Otherwise, it is set to 1. Defaults to FALSE.
#' @param fit_delay Should the delay of the plasma parent fraction be fitted? Otherwise, it is set to 0. Defaults to FALSE.
#' @param lower Named list of the lower limits.
#' @param upper Named list of the upper limits.
#' @param multstart_lower Named list of the lower starting limits.
#' @param multstart_upper Named list of the upper starting limits.
#' @param multstart_iter Number of fits to perform before deciding on an optimal.
#'
#' @return An nls fit object.
#' @export
#'
#' @examples
#' \dontrun{
#' pf <- bd_getdata(blooddata, output = "parentFraction")
#' metab_hill(pf$time, pf$parentFraction)
#' }
metab_hill <- function(time, parentFraction,
                       fit_ppf0 = FALSE,
                       fit_delay = FALSE,
                       lower = list(a = 0, b = 1, c = 0, ppf0 = 0.8, delay = -30),
                       upper = list(a = 1, b = 100, c = 100, ppf0 = 1.1, delay = 30),
                       multstart_lower = NULL,
                       multstart_upper = NULL,
                       multstart_iter = 100) {
  pf <- tibble::tibble(time = time, parentFraction = parentFraction)
  pf <- dplyr::arrange(pf, time)

  formula <- paste("parentFraction ~ metab_hill_model(time, a, b, c, ",
    "ppf0", ifelse(fit_ppf0, "", "=1"), ", ",
    "delay", ifelse(fit_delay, "", "=0"), ")",
    sep = ""
  )

  if (!fit_ppf0) {
    lower$ppf0 <- NULL
    upper$ppf0 <- NULL
  }

  if (!fit_delay) {
    lower$delay <- NULL
    upper$delay <- NULL
  }

  lower <- as.numeric(as.data.frame(lower))
  upper <- as.numeric(as.data.frame(upper))

  if (is.null(multstart_lower)) {
    multstart_lower <- lower
  } else {
    multstart_lower <- as.numeric(as.data.frame(multstart_lower))
  }

  if (is.null(multstart_upper)) {
    multstart_upper <- upper
  } else {
    multstart_upper <- as.numeric(as.data.frame(multstart_upper))
  }

  nls.multstart::nls_multstart(as.formula(formula),
    data = pf,
    lower = lower,
    upper = upper,
    start_lower = multstart_lower,
    start_upper = multstart_upper,
    iter = multstart_iter,
    supp_errors = "Y"
  )
}


# #' Extended Hill Function Model for Parent Fraction
# #'
# #' This is the model function for fitting of the extended Hill function.
# #'
# #' @param time Time in seconds.
# #' @param a Parameter A.
# #' @param b Parameter B.
# #' @param c Parameter C as the log10 of the 'classic' c parameter.
# #' @param d Parameter D.
# #' @param ppf0 The starting point of the parent fraction curve.
# #' @param delay The delay of the metabolism curve.
# #'
# #' @return Predicted values
# #' @export
# #'
# #' @examples
# #' metab_hillextended_model(seq(0, 60*60, by=120), 0.1, 3, 6.9, 0, 1, 0)
# metab_hillextended_model <- function(time, a, b, c, d, ppf0=1, delay=0) {
#
#   tcorr <- time-delay
#   t_before <- tcorr[ which(!(tcorr > 0)) ]
#   t_after <- tcorr[ which(tcorr > 0) ]
#
#   ind1_out <- rep(ppf0, length.out=length(t_before))
#   ind2_out <- ppf0 - (
#     ( ((ppf0-a) - d*t_after ) * t_after^b) /
#       ( 10^c + (t_after)^b )
#     )
#
#   out <- c(ind1_out, ind2_out)
#
#   return(out)
# }
#
#
# metab_hillextended <- function(time, parentFraction,
#                        fit_ppf0 = FALSE,
#                        fit_delay = FALSE,
#                        lower=list(a=0, b=1, c=0, d=-100, ppf0=0.8, delay=-30),
#                        upper=list(a=1, b=100, c=100, d=100, ppf0=1.1, delay=30),
#                        multstart_lower = NULL,
#                        multstart_upper = NULL,
#                        multstart_iter=100) {
#
#
#   pf <- tibble::tibble(time = time, parentFraction = parentFraction)
#   pf <- dplyr::arrange(pf, time)
#
#   formula <- paste("parentFraction ~ metab_hillextended_model(time, a, b, c, d, ",
#                    "ppf0", ifelse(fit_ppf0, "", "=1"), ", ",
#                    "delay", ifelse(fit_delay, "", "=0"), ")",
#                    sep = "")
#
#   if(!fit_ppf0) {
#     lower$ppf0 <- NULL
#     upper$ppf0 <- NULL
#   }
#
#   if(!fit_delay) {
#     lower$delay <- NULL
#     upper$delay <- NULL
#   }
#
#   lower <- as.numeric(as.data.frame(lower))
#   upper <- as.numeric(as.data.frame(upper))
#
#   if(is.null(multstart_lower)) {
#     multstart_lower <- lower
#   } else {
#   multstart_lower <- as.numeric(as.data.frame(multstart_lower))
#   }
#
#   if(is.null(multstart_upper)) {
#     multstart_upper <- upper
#   } else {
#     multstart_upper <- as.numeric(as.data.frame(multstart_upper))
#   }
#
#   nls.multstart::nls_multstart(as.formula(formula),
#                                data=pf,
#                                lower=lower,
#                                upper=upper,
#                                start_lower = multstart_lower,
#                                start_upper = multstart_upper,
#                                iter = multstart_iter,
#                                supp_errors = "Y")
#
# }



#' Hill Function Model modified by Guo et al for Parent Fraction
#'
#' This is the model function for fitting of the modification of the Hill function by Guo et al. (2013).
#'
#' @param time Time in seconds.
#' @param a Parameter A.
#' @param b Parameter B.
#' @param c Parameter C.
#' @param ppf0 The starting point of the parent fraction curve.
#' @param delay The delay of the metabolism curve.
#'
#' @return Predicted values
#' @export
#'
#' @references Guo Q, Colasanti A, Owen DR, et al. Quantification of the specific translocator protein signal of 18F-PBR111 in healthy humans: a genetic polymorphism effect on in vivo binding. J Nucl Med 2013; 54: 1915–1923.
#'
#' @examples
#' metab_hillguo_model(seq(0, 60 * 60, by = 120), 7, 0.6, 0.04, 1, 0)
metab_hillguo_model <- function(time, a, b, c, ppf0 = 1, delay = 0) {
  tcorr <- time - delay
  t_before <- tcorr[ which(!(tcorr > 0)) ]
  t_after <- tcorr[ which(tcorr > 0) ]

  ind1_out <- rep(ppf0, length.out = length(t_before))

  ind2_inner <- 1 - (t_after^3 / (t_after^3 + 10^a))
  ind2_out <- (ind2_inner^b + c) / (1 + c)

  out <- c(ind1_out, ind2_out)

  return(out)
}


#' Fit a Hill Function modification by Guo et al (2013) for Parent Fraction.
#'
#' This function fits the modification of the Hill function by Guo et al. (2013) to parent fraction data.
#'
#' @param time Time in seconds.
#' @param parentFraction Measured values of parent fraction.
#' @param fit_ppf0 Should the starting plasma parent fraction be fitted? Otherwise, it is set to 1. Defaults to FALSE.
#' @param fit_delay Should the delay of the plasma parent fraction be fitted? Otherwise, it is set to 0. Defaults to FALSE.
#' @param lower Named list of the lower limits.
#' @param upper Named list of the upper limits.
#' @param multstart_lower Named list of the lower starting limits.
#' @param multstart_upper Named list of the upper starting limits.
#' @param multstart_iter Number of fits to perform before deciding on an optimal.
#'
#' @return An nls fit object.
#' @export
#'
#' @references Guo Q, Colasanti A, Owen DR, et al. Quantification of the specific translocator protein signal of 18F-PBR111 in healthy humans: a genetic polymorphism effect on in vivo binding. J Nucl Med 2013; 54: 1915–1923.
#'
#' @examples
#' \dontrun{
#' pf <- bd_getdata(blooddata, output = "parentFraction")
#' metab_hillguo(pf$time, pf$parentFraction)
#' }
metab_hillguo <- function(time, parentFraction,
                          fit_ppf0 = FALSE,
                          fit_delay = FALSE,
                          lower = list(a = 0, b = 0, c = 0, ppf0 = 0.8, delay = -30),
                          upper = list(a = 100, b = 100, c = 100, ppf0 = 1.1, delay = 30),
                          multstart_lower = NULL,
                          multstart_upper = NULL,
                          multstart_iter = 100) {
  pf <- tibble::tibble(time = time, parentFraction = parentFraction)
  pf <- dplyr::arrange(pf, time)

  formula <- paste("parentFraction ~ metab_hillguo_model(time, a, b, c, ",
    "ppf0", ifelse(fit_ppf0, "", "=1"), ", ",
    "delay", ifelse(fit_delay, "", "=0"), ")",
    sep = ""
  )

  if (!fit_ppf0) {
    lower$ppf0 <- NULL
    upper$ppf0 <- NULL
  }

  if (!fit_delay) {
    lower$delay <- NULL
    upper$delay <- NULL
  }

  lower <- as.numeric(as.data.frame(lower))
  upper <- as.numeric(as.data.frame(upper))

  if (is.null(multstart_lower)) {
    multstart_lower <- lower
  } else {
    multstart_lower <- as.numeric(as.data.frame(multstart_lower))
  }

  if (is.null(multstart_upper)) {
    multstart_upper <- upper
  } else {
    multstart_upper <- as.numeric(as.data.frame(multstart_upper))
  }

  nls.multstart::nls_multstart(as.formula(formula),
    data = pf,
    lower = lower,
    upper = upper,
    start_lower = multstart_lower,
    start_upper = multstart_upper,
    iter = multstart_iter,
    supp_errors = "Y"
  )
}



#' Power Model for Parent Fraction
#'
#' This is the model function for fitting of the classic power function.
#'
#' @param time Time in seconds.
#' @param a Parameter A.
#' @param b Parameter B.
#' @param c Parameter C.
#' @param ppf0 The starting point of the parent fraction curve.
#' @param delay The delay of the metabolism curve.
#'
#' @return Predicted values
#' @export
#'
#' @examples
#' metab_power_model(seq(0, 60 * 60, by = 120), 0.004, 4.5, 0.27, 1, 0)
metab_power_model <- function(time, a, b, c, ppf0 = 1, delay = 0) {
  tcorr <- time - delay
  t_before <- tcorr[ which(!(tcorr > 0)) ]
  t_after <- tcorr[ which(tcorr > 0) ]

  ind1_out <- rep(ppf0, length.out = length(t_before))
  ind2_out <- ppf0 / (1 + (a * (t_after))^b)^c

  out <- c(ind1_out, ind2_out)

  return(out)
}


#' Fit a Power Function for Modelling Parent Fraction.
#'
#' This function fits the power function to parent fraction data.
#'
#' @param time Time in seconds.
#' @param parentFraction Measured values of parent fraction.
#' @param fit_ppf0 Should the starting plasma parent fraction be fitted? Otherwise, it is set to 1. Defaults to FALSE.
#' @param fit_delay Should the delay of the plasma parent fraction be fitted? Otherwise, it is set to 0. Defaults to FALSE.
#' @param lower Named list of the lower limits.
#' @param upper Named list of the upper limits.
#' @param multstart_lower Named list of the lower starting limits.
#' @param multstart_upper Named list of the upper starting limits.
#' @param multstart_iter Number of fits to perform before deciding on an optimal.
#'
#' @return An nls fit object.
#' @export
#'
#' @examples
#' \dontrun{
#' pf <- bd_getdata(blooddata, output = "parentFraction")
#' metab_power(pf$time, pf$parentFraction)
#' }
metab_power <- function(time, parentFraction,
                        fit_ppf0 = FALSE,
                        fit_delay = FALSE,
                        lower = list(a = 0, b = 1, c = 0, ppf0 = 0.8, delay = -30),
                        upper = list(a = 1, b = 10, c = 5, ppf0 = 1.1, delay = 30),
                        multstart_lower = NULL,
                        multstart_upper = NULL,
                        multstart_iter = 100) {
  pf <- tibble::tibble(time = time, parentFraction = parentFraction)
  pf <- dplyr::arrange(pf, time)

  formula <- paste("parentFraction ~ metab_power_model(time, a, b, c, ",
    "ppf0", ifelse(fit_ppf0, "", "=1"), ", ",
    "delay", ifelse(fit_delay, "", "=0"), ")",
    sep = ""
  )

  if (!fit_ppf0) {
    lower$ppf0 <- NULL
    upper$ppf0 <- NULL
  }

  if (!fit_delay) {
    lower$delay <- NULL
    upper$delay <- NULL
  }

  lower <- as.numeric(as.data.frame(lower))
  upper <- as.numeric(as.data.frame(upper))

  if (is.null(multstart_lower)) {
    multstart_lower <- lower
  } else {
    multstart_lower <- as.numeric(as.data.frame(multstart_lower))
  }

  if (is.null(multstart_upper)) {
    multstart_upper <- upper
  } else {
    multstart_upper <- as.numeric(as.data.frame(multstart_upper))
  }

  nls.multstart::nls_multstart(as.formula(formula),
    data = pf,
    lower = lower,
    upper = upper,
    start_lower = multstart_lower,
    start_upper = multstart_upper,
    iter = multstart_iter,
    supp_errors = "Y"
  )
}


#' Exponential Function Model for Parent Fraction
#'
#' This is the model function for fitting of the classic exponential function.
#'
#' @param time Time in seconds.
#' @param a Parameter A.
#' @param b Parameter B.
#' @param c Parameter C.
#' @param ppf0 The starting point of the parent fraction curve.
#' @param delay The delay of the metabolism curve.
#'
#' @return Predicted values
#' @export
#'
#' @examples
#' metab_exponential_model(seq(0, 60 * 60, by = 120), 0.02, 0, 0.001, 1, 0)
metab_exponential_model <- function(time, a, b, c, ppf0 = 1, delay = 0) {
  tcorr <- time - delay
  t_before <- tcorr[ which(!(tcorr > 0)) ]
  t_after <- tcorr[ which(tcorr > 0) ]

  ind1_out <- rep(ppf0, length.out = length(t_before))
  ind2_out <- a * exp(-b * t_after) + (ppf0 - a) * exp(-c * t_after)

  out <- c(ind1_out, ind2_out)

  return(out)
}



#' Fit an Exponential Function for Modelling Parent Fraction.
#'
#' This function fits the exponential function to parent fraction data.
#'
#' @param time Time in seconds.
#' @param parentFraction Measured values of parent fraction.
#' @param fit_ppf0 Should the starting plasma parent fraction be fitted? Otherwise, it is set to 1. Defaults to FALSE.
#' @param fit_delay Should the delay of the plasma parent fraction be fitted? Otherwise, it is set to 0. Defaults to FALSE.
#' @param lower Named list of the lower limits.
#' @param upper Named list of the upper limits.
#' @param multstart_lower Named list of the lower starting limits.
#' @param multstart_upper Named list of the upper starting limits.
#' @param multstart_iter Number of fits to perform before deciding on an optimal.
#'
#' @return An nls fit object.
#' @export
#'
#' @examples
#' \dontrun{
#' pf <- bd_getdata(blooddata, output = "parentFraction")
#' metab_exponential(pf$time, pf$parentFraction)
#' }
metab_exponential <- function(time, parentFraction,
                              fit_ppf0 = FALSE,
                              fit_delay = FALSE,
                              lower = list(a = 0, b = 0, c = 0, ppf0 = 0.8, delay = -30),
                              upper = list(a = 1, b = 1, c = 1, ppf0 = 1.1, delay = 30),
                              multstart_lower = NULL,
                              multstart_upper = NULL,
                              multstart_iter = 100) {
  pf <- tibble::tibble(time = time, parentFraction = parentFraction)
  pf <- dplyr::arrange(pf, time)

  formula <- paste("parentFraction ~ metab_exponential_model(time, a, b, c, ",
    "ppf0", ifelse(fit_ppf0, "", "=1"), ", ",
    "delay", ifelse(fit_delay, "", "=0"), ")",
    sep = ""
  )

  if (!fit_ppf0) {
    lower$ppf0 <- NULL
    upper$ppf0 <- NULL
  }

  if (!fit_delay) {
    lower$delay <- NULL
    upper$delay <- NULL
  }

  lower <- as.numeric(as.data.frame(lower))
  upper <- as.numeric(as.data.frame(upper))

  if (is.null(multstart_lower)) {
    multstart_lower <- lower
  } else {
    multstart_lower <- as.numeric(as.data.frame(multstart_lower))
  }

  if (is.null(multstart_upper)) {
    multstart_upper <- upper
  } else {
    multstart_upper <- as.numeric(as.data.frame(multstart_upper))
  }

  nls.multstart::nls_multstart(as.formula(formula),
    data = pf,
    lower = lower,
    upper = upper,
    start_lower = multstart_lower,
    start_upper = multstart_upper,
    iter = multstart_iter,
    supp_errors = "Y"
  )
}


#' Inverse Gamma Function Model for Parent Fraction
#'
#' This is the model function for fitting of the inverse Gamma function for the parent fraction.
#'
#' @param time Time in seconds.
#' @param shape Shape parameter.
#' @param rate Rate parameter.
#' @param ppf0 The starting point of the parent fraction curve.
#' @param delay The delay of the metabolism curve.
#'
#' @return Predicted values
#' @export
#'
#' @examples
#' metab_invgamma_model(seq(0, 60 * 60, by = 120), 1.97, 708, 1, 0)
metab_invgamma_model <- function(time, shape, rate, ppf0 = 1, delay = 0) {
  tcorr <- time - delay
  t_before <- tcorr[ which(!(tcorr > 0)) ]
  t_after <- tcorr[ which(tcorr > 0) ]

  ind1_out <- rep(ppf0, length.out = length(t_before))
  ind2_out <- ppf0 - invgamma::pinvgamma(t_after, shape, rate)

  out <- c(ind1_out, ind2_out)

  return(out)
}



#' Fit the Inverted Gamma Function for Modelling Parent Fraction.
#'
#' This function fits the inverted gamma function to parent fraction data.
#'
#' @param time Time in seconds.
#' @param parentFraction Measured values of parent fraction.
#' @param fit_ppf0 Should the starting plasma parent fraction be fitted? Otherwise, it is set to 1. Defaults to FALSE.
#' @param fit_delay Should the delay of the plasma parent fraction be fitted? Otherwise, it is set to 0. Defaults to FALSE.
#' @param lower Named list of the lower limits.
#' @param upper Named list of the upper limits.
#' @param multstart_lower Named list of the lower starting limits.
#' @param multstart_upper Named list of the upper starting limits.
#' @param multstart_iter Number of fits to perform before deciding on an optimal.
#'
#' @return An nls fit object.
#' @export
#'
#' @examples
#' \dontrun{
#' pf <- bd_getdata(blooddata, output = "parentFraction")
#' metab_invgamma(pf$time, pf$parentFraction)
#' }
metab_invgamma <- function(time, parentFraction,
                           fit_ppf0 = FALSE,
                           fit_delay = FALSE,
                           lower = list(shape = 0, rate = 0, ppf0 = 0.8, delay = -30),
                           upper = list(shape = 1000, rate = 1000, ppf0 = 1.1, delay = 30),
                           multstart_lower = NULL,
                           multstart_upper = NULL,
                           multstart_iter = 100) {
  pf <- tibble::tibble(time = time, parentFraction = parentFraction)
  pf <- dplyr::arrange(pf, time)

  formula <- paste("parentFraction ~ metab_invgamma_model(time, shape, rate, ",
    "ppf0", ifelse(fit_ppf0, "", "=1"), ", ",
    "delay", ifelse(fit_delay, "", "=0"), ")",
    sep = ""
  )

  if (!fit_ppf0) {
    lower$ppf0 <- NULL
    upper$ppf0 <- NULL
  }

  if (!fit_delay) {
    lower$delay <- NULL
    upper$delay <- NULL
  }

  lower <- as.numeric(as.data.frame(lower))
  upper <- as.numeric(as.data.frame(upper))

  if (is.null(multstart_lower)) {
    multstart_lower <- lower
  } else {
    multstart_lower <- as.numeric(as.data.frame(multstart_lower))
  }

  if (is.null(multstart_upper)) {
    multstart_upper <- upper
  } else {
    multstart_upper <- as.numeric(as.data.frame(multstart_upper))
  }

  nls.multstart::nls_multstart(as.formula(formula),
    data = pf,
    lower = lower,
    upper = upper,
    start_lower = multstart_lower,
    start_upper = multstart_upper,
    iter = multstart_iter,
    supp_errors = "Y"
  )
}
