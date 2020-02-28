#' Interpolate Blood Curves
#'
#' Function to interpolate the blood, plasma and parent fraction data in order
#' to create an \code{input} object for use with models requiring arterial input
#' functions. This function is a poor stand-in for the create_blooddata_*
#' functions.
#'
#' @param t_blood Numeric vector of times for each blood measurement in minutes.
#' @param blood Numeric vector of the radioactivity concentration in each blood
#'   measurement.
#' @param t_plasma Numeric vector of times for each plasma measurement in
#'   minutes.
#' @param plasma Numeric vector of the radioactivity concentration in each
#'   plasma measurement.
#' @param t_parentfrac Numeric vector of times for each parent fraction
#'   measurement in minutes.
#' @param parentfrac Numeric vector of the radioactivity concentration in each
#'   plasma measurement.
#' @param interpPoints The number of points to interpolate into.
#'
#'
#' @return A dataframe containing the time, blood, plasma and parent fraction
#'   interpolated into the same times, with \code{interpPoints} number of
#'   points.
#'
#' @description This function sorts out all the blood, plasma and parent
#'   fraction measurements into one convenient data frame for arterial models.
#'   It makes several 'editorial decisions' in the process. i) The data is
#'   interpolated into 6000 points by default (after some trial and error, the
#'   number of points can have dramatic implications for the speed of the
#'   function: 6000 and 1024 are very fast). ii) The different measurements are
#'   set to have the same times, thus if one measurement is taken for a shorter
#'   period than the others, the others are extended to that time point. This
#'   extension is performed by keeping the same value as the previous recorded
#'   value at that point. Again, get in touch if you want a better method
#'   introduced here. iii) This function sets the blood concentration, plasma
#'   concentration and parent fraction to 0, 0, and 1 respectively at time 0.
#'   Further, it removes any measurements at time <= 0. Best to add a bit to all
#'   measurements if you have time<=0 values, which can be fixed in the time
#'   shifting.
#'
#'
#' @examples
#' \dontrun{
#' input <- blood_interp(
#'   t_blood = blooddata$Time.sec. / 60,
#'   blood = blooddata$Cbl.nCi.cc.,
#'   t_plasma = plasmadata$Time.sec. / 60,
#'   plasma = plasmadata$Cpl.nCi.cc.,
#'   t_parentfrac = parentdata$Times / 60,
#'   parentfrac = parentdata$Fraction
#' )
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

blood_interp <- function(t_blood, blood, t_plasma, plasma, t_parentfrac, parentfrac, interpPoints = 6000) {
  if (max(t_blood) > 300 || max(t_plasma) > 300 || max(t_parentfrac) > 300) {
    warning("
            ***********************************************************************
            It looks like you have included seconds instead of minutes for time
            in at least one of the following: t_blood, t_plasma or t_parentfrac.
            This can cause wrong/weird results, and should be avoided. Ignore this
            warning if you just have really long measurements (over 300 minutes).
            ***********************************************************************")
  }

  blooddf <- data.frame(Time = t_blood, Value = blood, Measure = "Blood")
  plasmadf <- data.frame(Time = t_plasma, Value = plasma, Measure = "Plasma")
  parentfracdf <- data.frame(Time = t_parentfrac, Value = parentfrac, Measure = "ParentFrac")

  input <- rbind(blooddf, plasmadf, parentfracdf)

  maxtime <- max(input$Time)

  endpoints <- data.frame(
    Time = c(rep(maxtime, 3)),
    Value = c(tail(blooddf$Value, 1), tail(plasmadf$Value, 1), tail(parentfracdf$Value, 1)),
    Measure = c("Blood", "Plasma", "ParentFrac")
  )

  input <- subset(input, input$Time < maxtime)
  input <- rbind(input, endpoints)

  input <- subset(input, input$Time > 0)
  zeropoints <- data.frame(
    Time = rep(0, 3),
    Value = c(0, 0, 1),
    Measure = c("Blood", "Plasma", "ParentFrac")
  )
  input <- rbind(zeropoints, input)

  interptime <- pracma::linspace(0, maxtime, interpPoints)

  interpcurves <- plyr::dlply(
    input, "Measure",
    function(x) pracma::interp1(x$Time,
        x$Value,
        interptime,
        method = "linear"
      )
  )

  input <- tibble::tibble(
    Time = interptime,
    Blood = interpcurves$Blood,
    Plasma = interpcurves$Plasma,
    ParentFraction = interpcurves$ParentFrac,
    AIF = interpcurves$Plasma * interpcurves$ParentFrac
  )

  class(input) <- c("interpblood", class(input))

  return(input)
}


#' Shift timings of TAC and Input
#'
#' Function to shift backwards and forwards the timings of the TAC vector and input data frame
#' to make them consistent.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param inpshift The number of minutes by which the times of the input data frame should be adjusted.
#' @param shifttac Optional. Can the TAC be shifted. If the input is shifted negatively, it will then
#' contain negative time values. If \code{shifttac=TRUE}, the TAC will be made later, but if \code{shifttac=FALSE},
#' the blood, plasma and parent fraction with negative values will be removed. Default is TRUE.
#'
#' @return A list containing the following after time shifting: the times of the TAC \code{out$t_tac},
#' the TAC \code{out$tac}, the input dataframe \code{out$input}, the interpolated time \code{out$interptime}
#' (which is the same as the time in the input dataframe), and the interpolated TAC \code{out$i_tac}.
#'
#' @description This function sorts out all the time shifting of the input and TAC so that they are in the
#' same time. This function makes several 'editorial decisions'. i) If the TAC is shifted positively, all frames
#' will be shifted positively, but the time=0 frame will remain at time=0, i.e. there are no extra frames added.
#' Shifting of the TAC only occurs if the input is shifted negatively (and can be turned off using \code{shifttac=F}).
#' ii) If the input is shifted, and is subsequently shorter than the TAC, an extra measurement will be
#' added at \code{max(t_tac)} with the same value as the last measurement. iii) If the input is shifted positively,
#' all interpolated times will be shifted by the specified amount, but an extra measurement is added at time=0 of
#' 0,0,1 for blood, plasma and parent fraction respectively (followed by interpolation into default 6000 equally spaced time
#' intervals in the new time window) i.e. not the same process as for the TACs: I figure that since the blood changes
#' so quickly, this is likely more close to the true kinetics. Get in touch if you have suggestions for this.
#'
#'
#' @examples
#' \dontrun{
#' inpshift <- 0.25 # 15 seconds
#' newValues <- shift_timings(t_tac, tac, input, inpshift)
#' t_tac <- newValues$t_tac
#' tac <- newValues$tac
#' input <- newValues$input
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

shift_timings <- function(t_tac, tac, input, inpshift, shifttac = T) {
  interpPoints <- nrow(input)

  tacdf <- data.frame(Time = t_tac, Value = tac)

  if (min(tacdf$Time) > 0) {
    tacdf <- rbind(c(0, 0), tacdf)
  }

  if (min(input$Time) > 0) {
    input <- rbind(c(0, 0, 0, 1), input)
  }

  if (min(tacdf$Time) < 0) {
    stop("There are negative times in the TAC")
  }


  # Do the shift here!
  input$Time <- input$Time + inpshift

  if (min(input$Time) < 0) {
    if (shifttac == T) {
      beforezerodelay <- min(input$Time)
      tacdf$Time <- tacdf$Time - beforezerodelay
      tacdf[1, ] <- c(0, 0)
      input$Time <- input$Time - beforezerodelay
    } else {
      input <- subset(input, input$Time >= 0)
    }
  }

  if (min(input$Time) > 0) {
    # input$Time[1] <- 0
    input <- rbind(c(0, 0, 0, 1), input)
  }

  if (max(input$Time) > max(tacdf$Time)) {
    last_n <- which(input$Time > max(tacdf$Time))[1]
    input <- input[1:last_n, ]
    input$Time[last_n] <- max(tacdf$Time)
  } else if (max(input$Time) < max(tacdf$Time)) {
    input <- rbind(input, c(max(tacdf$Time), as.numeric(tail(input, 1)[-1])))
  }

  # Interpolating everything to the same time space

  interptime <- pracma::linspace(head(tacdf$Time, 1), tail(tacdf$Time, 1), interpPoints)

  input <- tibble::tibble(
    Time = interptime,
    Blood = pracma::interp1(input$Time, input$Blood, interptime,
      method = "linear"
    ),
    Plasma = pracma::interp1(input$Time, input$Plasma, interptime,
      method = "linear"
    ),
    ParentFraction = pracma::interp1(input$Time, input$ParentFraction, interptime,
      method = "linear"
    ),
    AIF = pracma::interp1(input$Time, input$AIF, interptime,
      method = "linear"
    )
  )

  class(input) <- c("interpblood", class(input))

  i_tac <- pracma::interp1(tacdf$Time, tacdf$Value, interptime, method = "linear")

  out <- list(t_tac = tacdf$Time, tac = tacdf$Value, input = input, interptime = interptime, i_tac = i_tac)
  return(out)
}


#' Shift timings of several TACs and Input
#'
#' Function to shift backwards and forwards the timings of the TAC vector and input data frame
#' to make them consistent. Similar to \code{shift_timings} but for a dataframe of TACs.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param tacsdf Dataframe of radioactivity concentrations in the target tissues for each frame in wide format, i.e. each
#' TAC should be a column. We include zero at time zero: if not included, it is added.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param inpshift The number of minutes by which the times of the input data frame should be adjusted.
#' @param shifttac Optional. Can the TACs be shifted. If the input is shifted negatively, it will then
#' contain negative time values. If \code{shifttac=TRUE}, the TACs will be made later, but if \code{shifttac=FALSE},
#' the blood, plasma and parent fraction with negative values will be removed. Default is TRUE.
#'
#' @return A list containing the following after time shifting: the times of the TACs \code{out$t_tac},
#' the TACs \code{out$tacdf}, the input dataframe \code{out$input}, the interpolated time \code{out$interptime}
#' (which is the same as the time in the input dataframe), and the interpolated TACs \code{out$i_tacdf}.
#'
#' @description This function sorts out all the time shifting of the input and TACs so that they are in the
#' same time. This function makes several 'editorial decisions'. i) If the TACs are shifted positively, all frames
#' will be shifted positively, but the time=0 frame will remain at time=0, i.e. there are no extra frames added.
#' Shifting of the TAC only occurs if the input is shifted negatively (and can be turned off using \code{shifttac=F}).
#' ii) If the input is shifted, and is subsequently shorter than the TAC, an extra measurement will be
#' added at \code{max(t_tac)} with the same value as the last measurement. iii) If the input is shifted positively,
#' all interpolated times will be shifted by the specified amount, but an extra measurement is added at time=0 of
#' 0,0,1 for blood, plasma and parent fraction respectively (followed by interpolation into 6000 equally spaced time
#' intervals in the new time window) i.e. not the same process as for the TACs: I figure that since the blood changes
#' so quickly, this is likely more close to the true kinetics. Get in touch if you have suggestions for this.
#'
#'
#' @examples
#' \dontrun{
#' inpshift <- 0.25 # 15 seconds
#' newValues <- shift_timings(t_tac, tac, input, inpshift)
#' t_tac <- newValues$t_tac
#' tac <- newValues$tac
#' input <- newValues$input
#' }
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

shift_timings_df <- function(t_tac, tacsdf, input, inpshift, shifttac = T) {
  interpPoints <- nrow(input)

  tacdf <- data.frame(Time = t_tac)
  tacdf <- cbind(tacdf, tacsdf)


  if (min(tacdf$Time) > 0) {
    tacdf <- rbind(rep(0, ncol(tacdf)), tacdf)
  }

  input$Time <- input$Time + inpshift

  if (min(input$Time) < 0) {
    if (shifttac == T) {
      beforezerodelay <- min(input$Time)
      tacdf$Time <- tacdf$Time - beforezerodelay
      tacdf[1, ] <- 0
      input$Time <- input$Time - beforezerodelay
    } else {
      input <- subset(input, input$Time >= 0)
    }
  }

  if (min(input$Time) > 0) {
    input <- rbind(c(0, 0, 0, 1), input)
  }

  if (max(input$Time) > max(tacdf$Time)) {
    last_n <- which(input$Time > max(tacdf$Time))[1]
    input <- input[1:last_n, ]
    input$Time[last_n] <- max(tacdf$Time)
  } else if (max(input$Time) < max(tacdf$Time)) {
    input <- rbind(input, c(max(tacdf$Time), as.numeric(tail(input, 1)[2:4])))
  }


  # Interpolating everything to the same time space

  interptime <- pracma::linspace(head(tacdf$Time, 1), tail(tacdf$Time, 1), interpPoints)


  input <- tibble::tibble(
    Time = interptime,
    Blood = pracma::interp1(input$Time, input$Blood, interptime,
      method = "linear"
    ),
    Plasma = pracma::interp1(input$Time, input$Plasma, interptime,
      method = "linear"
    ),
    ParentFraction = pracma::interp1(input$Time, input$ParentFraction, interptime,
      method = "linear"
    ),
    AIF = pracma::interp1(input$Time, input$AIF, interptime,
      method = "linear"
    )
  )

  class(input) <- c("interpblood", class(input))

  i_tacs <- sapply(tacsdf, function(x) pracma::interp1(tacdf$Time, x, interptime, method = "linear"))

  out <- list(
    t_tac = tacdf$Time, tacdf = tacdf[, -1], input = input, interptime = interptime,
    i_tacdf = i_tacs
  )
  return(out)
}


#' Plot the Timings of the TAC and Arterial Input Function
#'
#' Function to compare the timings of the the TAC and the Arterial Input Function (AIF) to see whether they are aligned.
#' Can be used for assessing the fit of the delay (\code{inpshift}).
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a
#' zero. If a time zero frame is not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time
#' zero: if not included, it is added.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param inpshift The number of minutes by which the times of the input data frame should be adjusted.
#' @param zoomTime The number of minutes to show for the close-up of the match between the TAC and AIF.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' data(pbr28)
#'
#' t_tac <- pbr28$tacs[[2]]$Times / 60
#' tac <- pbr28$tacs[[2]]$FC
#'
#' input <- blood_interp(
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cbl_dispcorr,
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cpl_metabcorr,
#'   t_parentfrac = 1, parentfrac = 1
#' )
#'
#' plot_inptac_timings(t_tac, tac, input, inpshift = 0.12, zoomTime = 5)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_inptac_timings <- function(t_tac, tac, input, inpshift, zoomTime = 5) {
  newvals <- shift_timings(t_tac, tac, input, inpshift)

  inp <- newvals$input$Plasma * newvals$input$ParentFraction
  interptime <- newvals$input$Time
  i_tac <- newvals$i_tac

  i_df <- data.frame(interptime, TAC = i_tac, AIF = inp)

  i_df_tidy <- tidyr::gather(
    i_df,
    key = Curve, value = Radioactivity,
    TAC, AIF
  )

  Curve <- factor(i_df_tidy$Curve)

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(Curve)
  colScale <- scale_colour_manual(name = "Curve", values = myColors)

  outplot <- ggplot(i_df_tidy, aes(x = interptime, y = Radioactivity, colour = Curve)) +
    geom_line() + xlab("Time (min)") +
    colScale + coord_cartesian(ylim = c(0, max(tac) * 1.5), xlim = c(0, zoomTime))

  return(outplot)
}



#' Plot the Timings of the TAC and Arterial Input Function from a Fit Object
#'
#' Function to compare the timings of the the TAC and the Arterial Input Function (AIF) to see whether they are aligned
#' from a kinetic fit object.  The fit needs to be a model utilising AIF (i.e. cannot be a reference region)
#'
#' @param fitout The output object of a fitting procedure which utilises arterial input function, e.g. \code{twotcm}.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#' @param zoomTime Optional. The number of minutes to show for the match between the TAC and AIF. Default is 5.
#'
#' @return A ggplot2 object of the plot.
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
#' twotcmout <- twotcm(t_tac, tac, input, weights, frameStartEnd = c(1, 25))
#'
#' plot_inptac_fit(twotcmout)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export

plot_inptac_fit <- function(fitout, roiname = NULL, zoomTime = 5) {
  if (is.null(roiname)) {
    roiname <- "ROI"
  }

  if (!("input" %in% names(fitout))) {
    stop("No input object found in fit output")
  }

  measureddf <- data.frame(
    Time = fitout$tacs$Time,
    Radioactivity = fitout$tacs$Target,
    Region = paste0(roiname, ".Measured")
  )

  inputdf <- data.frame(
    Time = fitout$input$Time,
    Radioactivity = fitout$input$Plasma * fitout$input$ParentFraction,
    Region = "AIF"
  )

  plotdf <- rbind(inputdf, measureddf)

  plotdf$Region <- forcats::fct_inorder(factor(plotdf$Region))

  myColors <- RColorBrewer::brewer.pal(3, "Set1")
  names(myColors) <- levels(plotdf$Region)
  colScale <- scale_colour_manual(name = "Region", values = myColors)

  outplot <- ggplot(plotdf, aes(x = Time, y = Radioactivity, colour = Region)) + colScale +
    geom_point(data = subset(plotdf, plotdf$Region == paste0(roiname, ".Measured")), aes(shape = "a")) +
    geom_line() +
    guides(shape = FALSE, color = guide_legend(order = 1)) +
    coord_cartesian(ylim = c(0, max(measureddf$Radioactivity) * 1.5), xlim = c(0, zoomTime))

  return(outplot)
}


#' Plot input
#'
#' Provides a plot of the input data to get an idea of how it looks. Blood, plasma and AIF data are scaled to a common value.
#'
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#'
#' @return A ggplot2 object
#'
#' @examples
#'
#' data(pbr28)
#'
#' input <- blood_interp(
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cbl_dispcorr,
#'   pbr28$procblood[[2]]$Time / 60, pbr28$procblood[[2]]$Cpl_metabcorr,
#'   t_parentfrac = 1, parentfrac = 1
#' )
#'
#' plot_input(input)
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @import ggplot2
#'
#' @export
plot_input <- function(input) {
  input$BPR <- input$Blood / input$Plasma
  input$Blood <- input$Blood / max(c(input$Blood, input$Plasma))
  input$Plasma <- input$Plasma / max(c(input$Blood, input$Plasma))
  input$AIF <- input$AIF

  input <- dplyr::rename(input,
    "Parent Fraction" = ParentFraction,
    "Blood-Plasma Ratio" = BPR
  )

  tidyinput <- tidyr::gather(input, Data, Radioactivity, -Time) %>%
    dplyr::mutate(Data = forcats::as_factor(Data))

  ggplot(tidyinput, aes(x = Time, y = Radioactivity, colour = Data)) +
    geom_line() +
    coord_cartesian(ylim = c(0, 1))
}

#' Perform dispersion correction on blood data collected using ABSS
#'
#' Perform dispersion correction on continuous blood samples from vectors as per
#' http://www.turkupetcentre.net/petanalysis/input_dispersion.html. The measured
#' values are first interpolated so that all measurements have an equal time
#' between them, so any missing measurements are interpolated.
#'
#' @param time Vector of time (in seconds) of each sample.
#' @param activity Vector of measured radioactivity concentrations in  blood
#'   samples.
#' @param tau Time constant denoting the time of dispersion (in seconds).
#' @param timedelta The time difference between each measured sample. Defaults
#'   to the most common time difference between the first 20 measurements.
#' @param keep_interpolated Defaults to TRUE: should interpolated samples which
#'   were not included in the original input be included in the output.
#' @param smooth_iterations The number of times that the smoothing of each pair
#'   of observations should be performed using blood_smooth. Defaults to 0.
#'
#' @return A tibble containing the measured times and radioactivity
#'   concentrations after dispersion correction.
#' @export
#'
#' @examples
#' time <- 1:20
#' activity <- rnorm(20)
#' tau <- 2.5
#'
#' blood_dispcor(time, activity, tau, 1)
blood_dispcor <- function(time, activity, tau, timedelta = NULL,
                          keep_interpolated = T, smooth_iterations = 0) {
  if (is.null(timedelta)) {
    diffs <- diff(head(time, 20))
    deltas <- sort(table(diffs), decreasing = TRUE)
    timedelta <- as.numeric(names(deltas)[1])
  }

  # linear interpolation of time and blood into seconds, in case values missing
  i_time <- seq(from = min(time), to = max(time), by = timedelta)
  i_activity <- pracma::interp1(time, activity, i_time, method = "linear")

  # Adding a value on each side to make the output the same length
  i_time <- c(head(i_time, 1), i_time, tail(i_time, 1))
  i_activity <- c(head(i_activity, 1) - timedelta, i_activity, tail(i_activity, 1) + timedelta)

  i_integ_activity <- pracma::cumtrapz(i_time, i_activity)
  i_integ_true <- tau * i_activity + i_integ_activity

  ind <- 2:(length(i_time) - 1)

  time_out <- i_time[ind]
  activity_out <- (i_integ_true[ind + 1] - i_integ_true[ind - 1]) /
    (2 * timedelta)


  if (keep_interpolated) {

    out <- tibble::tibble(time = time_out,
                          activity = activity_out)
  } else {

    activity_out <- pracma::interp1(x = time_out,
                               y = activity_out,
                               xi = time)

    out <- tibble::tibble(time = time,
                          activity = activity_out)
  }

  out <- blood_smooth(out$time, out$activity, smooth_iterations)

  return(out)
}



#' Smooth continuous blood data
#'
#' This function averages each pair of blood measurements from ABSS systems -
#' can be used after dispersion correction.
#'
#' @param time Vector of time (in seconds) of each sample.
#' @param activity Vector of measured radioactivity concentrations in  blood
#'   samples.
#' @param iterations The number of times that the smoothing of each pair of
#'   observations should be performed.
#'
#' @return A tibble containing the measured times and radioactivity
#'   concentrations after smoothing.
#' @export
#'
#' @examples
#' time <- 1:20
#' activity <- rnorm(20)
#'
#' blood_smooth(time, activity, 5)
blood_smooth <- function(time, activity, iterations = 1) {
  if (iterations > 0) {
    for (i in 1:iterations) {
      idx <- 1:(length(time) - 1)

      time <- (time[idx] + time[idx + 1]) / 2
      activity <- (activity[idx] + activity[idx + 1]) / 2
    }
  }

  tibble::tibble(time = time, activity = activity)
}
