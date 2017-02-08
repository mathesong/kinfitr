#' Interpolate Blood Curves
#'
#' Function to interpolate the blood, plasma and parent fraction data in order to create an 
#' \code{input} object for use with models requiring arterial input functions.
#'
#' @param t_blood Numeric vector of times for each blood measurement in minutes.
#' @param blood Numeric vector of the radioactivity concentration in each blood measurement.
#' @param t_plasma Numeric vector of times for each plasma measurement in minutes.
#' @param plasma Numeric vector of the radioactivity concentration in each plasma measurement.
#' @param t_parentfrac Numeric vector of times for each parent fraction measurement in minutes.
#' @param parentfrac Numeric vector of the radioactivity concentration in each plasma measurement.
#'
#'
#' @return A dataframe containing the time, blood, plasma and parent fraction interpolated into the 
#' same times, with \code{interpPoints} number of points.
#'
#' @description This function sorts out all the blood, plasma and parent fraction measurements into
#' one convenient data frame for arterial models. It makes several 'editorial decisions' in the process. 
#' i) The data is interpolated into 4096 points. This is a good number of points for Fast Fourier Transform as
#' it is a power of 2, and it is sufficiently many to have a sufficient level of detail. Get in touch if there
#' is a good reason to change this value. ii) The different measurements are set to have the same times, thus
#' if one measurement is taken for a shorter period than the others, the others are extended to that time point.
#' This extension is performed by keeping the same value as the previous recorded value at that point. Again, get
#' in touch if you want a better method introduced here. iii) This function sets the blood concentration, plasma
#' concentration and parent fraction to 0, 0, and 1 respectively at time 0. Further, it removes any measurements
#' at time <= 0. Best to add a bit to all measurements if you have time<=0 values, which can be fixed in the time
#' shifting.
#' 
#'
#' @examples
#' input <- blood_interp(t_blood = blooddata$Time.sec./60, 
#'            blood = blooddata$Cbl.nCi.cc., 
#'            t_plasma = plasmadata$Time.sec./60, 
#'            plasma = plasmadata$Cpl.nCi.cc., 
#'            t_parentfrac = parentdata$Times/60, 
#'            parentfrac = parentdata$Fraction)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

blood_interp <- function(t_blood, blood, t_plasma, plasma, t_parentfrac, parentfrac) {
  
  if( max(t_blood)>300 || max(t_plasma)>300 || max(t_parentfrac)>300 ) {
    warning('
      ******************************************************************************\n
          It looks like you have included seconds instead of minutes for time\n
          in at least one of the following: t_blood, t_plasma or t_parentfrac.\n
          This can cause wrong/weird results, and should be avoided. Ignore this\n
          warning if you just have really long measurements (over 300 minutes).\n
      ******************************************************************************')
  }
  
  interpPoints = 4096
  
  blooddf <- data.frame(Time=t_blood, Value=blood, Measure='Blood')
  plasmadf <- data.frame(Time=t_plasma, Value=plasma, Measure='Plasma')
  parentfracdf <- data.frame(Time=t_parentfrac, Value=parentfrac, Measure='ParentFrac')
  
  input <- rbind(blooddf, plasmadf, parentfracdf)
  
  maxtime = max(input$Time)
  
  endpoints <- data.frame(Time = c( rep(maxtime, 3) ) , 
                          Value = c( tail(blooddf$Value, 1) , tail(plasmadf$Value, 1) , tail(parentfracdf$Value, 1) ) ,
                          Measure = c( 'Blood', 'Plasma', 'ParentFrac' )    )
  
  input <- subset(input, input$Time < maxtime)
  input <- rbind(input, endpoints)
  
  input <- subset(input, input$Time > 0)
  zeropoints <- data.frame(Time = rep(0, 3) , 
                           Value = c(0, 0, 1) ,
                           Measure = c( 'Blood', 'Plasma','ParentFrac' )    )
  input <- rbind(zeropoints, input)
  
  interptime <- pracma::linspace( 0 , maxtime , interpPoints )
  
  interpcurves <- plyr::dlply(input, 'Measure', function(x) pracma::interp1(x$Time, x$Value, interptime, method="linear"))
  
  interped <- data.frame(Time=interptime, blood=interpcurves$Blood, plasma=interpcurves$Plasma, parentfrac=interpcurves$ParentFrac)
  
  return(interped)
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
#' 0,0,1 for blood, plasma and parent fraction respectively (followed by interpolation into 4096 equally spaced time 
#' intervals in the new time window) i.e. not the same process as for the TACs: I figure that since the blood changes
#' so quickly, this is likely more close to the true kinetics. Get in touch if you have suggestions for this.
#' 
#'
#' @examples
#' inpshift = 0.25   # 15 seconds
#' newValues <- shift_timings(t_tac, tac, input, inpshift)
#' t_tac <- newValues$t_tac
#' tac <- newValues$tac
#' input <- newValues$input
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

shift_timings <- function(t_tac, tac, input, inpshift, shifttac=T) {
  
  interpPoints = 4096
  
  tacdf <- data.frame(Time = t_tac, Value = tac)
  
  if(min(tacdf$Time) > 0) {
    tacdf <- rbind(c(0,0), tacdf)
  }
  
  if(min(input$Time) > 0) {
    input <- rbind(c(0,0,0,1), input)
  }
  
  if(min(tacdf$Time) < 0) {
    stop('There are negative times in the TAC')
  }
  
  
  # Do the shift here!
  input$Time <- input$Time + inpshift
  
  if(min(input$Time) < 0) {
    
    if(shifttac == T) {
      beforezerodelay = min(input$Time)
      tacdf$Time <- tacdf$Time-beforezerodelay
      tacdf[1,] = c(0,0)
      input$Time <- input$Time-beforezerodelay
    } else {
      input <- subset(input, input$Time >=0)
    }
  }
  
  if(min(input$Time) > 0) {
    #input$Time[1] <- 0
    input <- rbind(c(0,0,0,1), input)
  }
  
  if(max(input$Time) > max(tacdf$Time)) {
    last_n <- which( input$Time > max(tacdf$Time) )[1]
    input <- input[1:last_n,]
    input$Time[last_n] <- max(tacdf$Time)
  } else if(max(input$Time) < max(tacdf$Time)) {
    input <- rbind( input, c(max(tacdf$Time),  as.numeric(tail(input, 1)[2:4])))
  }
  
  # Interpolating everything to the same time space
  
  interptime <- pracma::linspace( head(tacdf$Time,1) , tail(tacdf$Time,1) , interpPoints )
  
  i_blood <- pracma::interp1(input$Time, input$blood, interptime, method="linear")
  i_plasma <- pracma::interp1(input$Time, input$plasma, interptime, method="linear")
  i_parentfrac <- pracma::interp1(input$Time, input$parentfrac, interptime, method="linear")
  i_tac <- pracma::interp1(tacdf$Time, tacdf$Value, interptime, method="linear")
  
  input <- data.frame(Time=interptime, blood=i_blood, plasma=i_plasma, parentfrac=i_parentfrac)
  
  
  
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
#' 0,0,1 for blood, plasma and parent fraction respectively (followed by interpolation into 4096 equally spaced time 
#' intervals in the new time window) i.e. not the same process as for the TACs: I figure that since the blood changes
#' so quickly, this is likely more close to the true kinetics. Get in touch if you have suggestions for this.
#' 
#'
#' @examples
#' inpshift = 0.25   # 15 seconds
#' newValues <- shift_timings(t_tac, tac, input, inpshift)
#' t_tac <- newValues$t_tac
#' tac <- newValues$tac
#' input <- newValues$input
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#'
#' @export

shift_timings_df <- function(t_tac, tacsdf, input, inpshift, shifttac=T) {
  
  interpPoints = 4096
  
  tacdf <- data.frame(Time = t_tac)
  tacdf <- cbind(tacdf, tacsdf)
  
  input$Time <- input$Time + inpshift
  
  if(min(input$Time) < 0) {
    
    if(shifttac == T) {
      beforezerodelay = min(input$Time)
      tacdf$Time <- tacdf$Time-beforezerodelay
      tacdf[1,] = 0
      input$Time <- input$Time-beforezerodelay
    } else {
      input <- subset(input, input$Time >=0)
    }
  }
  
  if(min(input$Time) > 0) {
    input <- rbind(c(0,0,0,1), input)
  }
  
  if(max(input$Time) > max(tacdf$Time) ) {
    last_n <- which( input$Time > max(tacdf$Time) )[1]
    input <- input[1:last_n,]
    input$Time[last_n] <- max(tacdf$Time)
  } else if(max(input$Time) < max(tacdf$Time)) {
    input <- rbind( input, c(max(tacdf$Time),  as.numeric(tail(input, 1)[2:4])))
  }
  
  
  # Interpolating everything to the same time space
  
  interptime <- pracma::linspace( head(tacdf$Time,1) , tail(tacdf$Time,1) , interpPoints )
  
  i_blood <- pracma::interp1(input$Time, input$blood, interptime, method="linear")
  i_plasma <- pracma::interp1(input$Time, input$plasma, interptime, method="linear")
  i_parentfrac <- pracma::interp1(input$Time, input$parentfrac, interptime, method="linear")
  i_tacs <- sapply(tacsdf, function(x) pracma::interp1(tacdf$Time, x, interptime, method="linear"))
  
  input <- data.frame(Time=interptime, blood=i_blood, plasma=i_plasma, parentfrac=i_parentfrac)
  
  out <- list(t_tac = tacdf$Time, tacdf = tacdf[,-1], input = input, interptime = interptime,
              i_tacdf = i_tacs)
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
#' @return Plots the two-panel figure 
#'
#' @examples
#' plot_inptac_timings(t_tac, tac, input, inpshift=0.25, zoomTime=5)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#' 
#' @import ggplot2
#'
#' @export

plot_inptac_timings <- function(t_tac, tac, input, inpshift, zoomTime=5) {
  
  newvals <- shift_timings(t_tac, tac, input, inpshift)
  
  # t_tac <- newvals$t_tac
  # tac <- newvals$tac
  inp <- newvals$input$plasma * newvals$input$parentfrac
  interptime <- newvals$input$Time
  i_tac <- newvals$i_tac
  
  # if(min(t_tac) < 0) {
  #   stop('There are negative times in the TAC')
  # }
  # 
  # if(min(t_tac) > 0) {
  #   t_tac = c(0, t_tac)
  #   tac = c(0, tac)
  # }
  # 
  # tacdf <- data.frame(Time = t_tac, Value = tac)
  # input <- data.frame(Time = t_inp, Value = inp)
  # 
  # input$Time <- input$Time + inpshift
  # 
  # if(min(input$Time) < 0) {
  #   beforezerodelay = min(input$Time)
  #   
  #   tacdf$Time <- tacdf$Time-beforezerodelay
  #   tacdf <- rbind(c(0,0), tacdf)
  #   
  #   input$Time <- input$Time-beforezerodelay
  # }
  
  # if(min(input$Time) > 0) {
  #   input <- rbind(c(0,0), input)
  # }
  # 
  # input <- subset(input, input$Time < max(tacdf$Time))
  # input <- rbind(input, c(max(tacdf$Time), tail(input$Value,1)))
  
  # i_time <- seq(0, max(tacdf$Time), by=max(tacdf$Time)/4096)
  # i_tac <- pracma::interp1(tacdf$Time, tacdf$Value, i_time)
  # i_inp <- pracma::interp1(input$Time, input$Value, i_time)
  
  # i_tac <- pracma::interp1(t_tac, tac, interptime)
  
  i_df <- data.frame(interptime, TAC = i_tac, AIF = inp)
  
  i_df_tidy <- tidyr::gather(i_df, key=Curve, value=Radioactivity, 
                             TAC, AIF)
  
  Curve <- factor(i_df_tidy$Curve)
  
  myColors <- RColorBrewer::brewer.pal(3,"Set1")
  names(myColors) <- levels(Curve)
  colScale <- scale_colour_manual(name = "Curve",values = myColors)
  
  inptiming_overall <- ggplot(i_df_tidy, aes(x = interptime, y = Radioactivity, colour = Curve)) + 
    geom_line() + xlim(0,max(t_tac)) + ylim(0,max(inp)) + xlab('Time (min)') +
    colScale
  inptiming_close <- ggplot(i_df_tidy, aes(x = interptime, y = Radioactivity, colour = Curve)) + 
    geom_line() + ylim(0,max(inp)/2) + xlab('Time (min)') +
    scale_x_continuous(breaks=seq(0,zoomTime, by = 1), limits = c(0,zoomTime)) +
    colScale
  
  gridExtra::grid.arrange(inptiming_overall, inptiming_close, ncol=1)
  
  inptiming <- gridExtra::arrangeGrob(inptiming_overall, inptiming_close, ncol=1)
  #return(inptiming)
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
#' plot_inptac_fit(twotcmout)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#' 
#' @import ggplot2
#'
#' @export

plot_inptac_fit <- function(fitout, roiname, zoomTime=5) {
  
  if(missing(roiname)) {roiname = 'ROI'}
  
  if(!('input' %in% names(fitout))) {
    stop('No input object found in fit output')
  }
  
  measureddf <- data.frame(Time = fitout$tacs$Time,
                           Radioactivity = fitout$tacs$Target,
                           Region=paste0(roiname, '.Measured'))
  
  inputdf <- data.frame(Time = fitout$input$Time,
                        Radioactivity = fitout$input$plasma*fitout$input$parentfrac,
                        Region='AIF')
  
  plotdf <- rbind(inputdf, measureddf)
  plotdf <- dplyr::filter(plotdf, Time < zoomTime)
  
  plotdf$Region <- forcats::fct_inorder(factor(plotdf$Region) )
  
  myColors <- RColorBrewer::brewer.pal(3,"Set1")
  names(myColors) <- levels( plotdf$Region )
  colScale <- scale_colour_manual(name = "Region",values = myColors)
  
  outplot = ggplot(plotdf, aes(x=Time, y=Radioactivity, colour=Region)) + colScale + 
    geom_point(data=subset(plotdf, plotdf$Region == paste0(roiname, '.Measured')), aes(shape='a')) + 
    geom_line() + 
    guides(shape=FALSE, color=guide_legend(order=1)) + scale_size(range=c(1,3)) + ylim(c(0,max(measureddf$Radioactivity)*1.5))
  
  return(outplot)
}