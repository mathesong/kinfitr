#' Patlak Reference Tissue Model
#'
#' Function to fit the Patlak Reference Tissue Model of Patlak & Blasbert (1985) to data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a 
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame. We include zero at 
#' time zero: if not included, it is added.
#' @param roitac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time 
#' zero: if not included, it is added.
#' @param tstarIncludedFrames The number of frames to be used in the regression model, i.e. the number of frames for which 
#' the function is linear after pseudo-equilibrium is reached. This is a count from the end of the measurement, so a value of 
#' 10 means that last 10 frames will be used. This value can be estimated using \code{refPatlak_tstar}.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero: 
#' if not included, it is added. If not specified, uniform weights will be used.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20). 
#' This is to assess time stability.
#'
#' @return A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit}, a dataframe 
#' containing the TACs of the data \code{out$tacs}, a dataframe containing the TACs of the fitted values \code{out$fitvals}, 
#' a vector of the weights \code{out$weights}, and the specified tstarIncludedFrames value \code{out$tstarIncludedFrames}
#'
#' @examples
#' refPatlak(t_tac, reftac, roitac, tstarIncludedFrames=10, weights=weights)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#' 
#' @references Patlak CS, Blasberg RG. Graphical evaluation of blood-to-brain transfer constants from multiple-time uptake data. Generalizations. Journal of Cerebral Blood Flow & Metabolism. 1985 Dec 1;5(4):584-90.
#'
#' @export

refPatlak <- function(t_tac, reftac, roitac, tstarIncludedFrames, weights, frameStartEnd) {
  
  
  # Tidying
  
  tidyinput <- tidyinput_ref(t_tac, reftac, roitac, weights, frameStartEnd)
  
  t_tac   <- tidyinput$t_tac
  reftac  <- tidyinput$reftac
  roitac  <- tidyinput$roitac
  weights <- tidyinput$weights
  
  
  # Parameters
  
  patlak_roi = roitac / reftac
  patlak_ref = pracma::cumtrapz(t_tac, reftac) / reftac
  
  patlak_equil_roi = tail(patlak_roi, tstarIncludedFrames)
  patlak_equil_ref = tail(patlak_ref, tstarIncludedFrames)
  weights_equil   = tail(weights,   tstarIncludedFrames)
  
  
  # Solution
  
  patlak_model <- lm( patlak_equil_roi ~ patlak_equil_ref, weights = weights_equil )
  
  # Output
  
  par = as.data.frame(list(K = as.numeric(patlak_model$coefficients[2])))
  fit = patlak_model
  
  tacs = data.frame(Time = t_tac, Reference = reftac, Target = roitac)
  
  fitvals = data.frame(Patlak_ROI = patlak_roi, Patlak_Ref = patlak_ref)
  
  out <- list( par = par, fit = fit, tacs = tacs, fitvals = fitvals, weights=weights, 
               tstarIncludedFrames = tstarIncludedFrames, model='refPatlak')
  
  return(out)
  
}

#' Plot: Patlak Reference Tissue Model
#'
#' Function to visualise the fit of the refPatlak model to data.
#'
#' @param refpatlakout The output object of the refPatlak fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' plot_refPatlakfit(refpatlakout)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#' 
#' @import ggplot2
#'
#' @export

plot_refPatlakfit <- function(refpatlakout, roiname = NULL) {
 
  plotdf <- data.frame(Weights = refpatlakout$weights,
                       Patlak_ref = refpatlakout$fitvals$Patlak_Ref,
                       Patlak_roi= refpatlakout$fitvals$Patlak_ROI,
                       Equilibrium = as.character('Before'))

  
  plotdf$Equilibrium <- as.character(plotdf$Equilibrium)
  plotdf$Equilibrium [ (nrow(plotdf)-(refpatlakout$tstarIncludedFrames-1)) :nrow(plotdf)  ] = 'After'
  
  plotdf$Equilibrium <- forcats::fct_inorder(factor(plotdf$Equilibrium))
  
  myColors <- RColorBrewer::brewer.pal(3,"Set1")
  names(myColors) <- levels(plotdf$Equilibrium)
  colScale <- scale_colour_manual(name = paste0(roiname, "\nEquilibrium"),values = myColors)
  
  outplot <- ggplot(data=plotdf, aes(x=Patlak_ref, y = Patlak_roi, colour=Equilibrium)) + 
    geom_point(aes(shape='a', size=Weights)) +
    geom_abline(slope     = as.numeric(refpatlakout$fit$coefficients[2]), 
                intercept = as.numeric(refpatlakout$fit$coefficients[1])) +
    xlab('Integ(C_Ref) / C_Ref') + ylab('C_Tissue / C_Ref') + colScale +
    guides(shape=FALSE, color=guide_legend(order=1)) + scale_size(range=c(1,3))
  
  return(outplot)
   
}


#' Tstar Finder: Patlak Reference Tissue Model
#'
#' Function to identify where t* is for the Patlak Reference Tissue Model.
#'
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a 
#' zero. If a time zero frame is not included, it will be added.
#' @param reftac Numeric vector of radioactivity concentrations in the reference tissue for each frame.
#' @param lowroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with low binding.
#' @param medroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with medium binding.
#' @param highroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with high binding.
#' @param filename The name of the output image: filename_refPatlak.jpeg
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_refPatlak.jpeg
#'
#' @examples
#' refPatlak_tstar(t_tac, reftac, taclow, tacmed, tachigh, 'demonstration')
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#' 
#' @import ggplot2
#'
#' @export

refPatlak_tstar <- function(t_tac, reftac, lowroi, medroi, highroi, filename, gridbreaks = 2) {
  
  frames = length(reftac)
  lowroi_fit  <- refPatlak(t_tac, reftac, lowroi, length(reftac))
  medroi_fit  <- refPatlak(t_tac, reftac, medroi, length(reftac))
  highroi_fit <- refPatlak(t_tac, reftac, highroi, length(reftac))
  
  patlak_xlab = 'Integ(C_Ref) / C_Ref'
  patlak_ylab = 'C_Tissue / C_Ref'
  
  low_linplot  <- qplot(lowroi_fit$fitvals$Patlak_Ref, lowroi_fit$fitvals$Patlak_ROI) + ggtitle('Low') + xlab(patlak_xlab) + ylab(patlak_ylab)
  med_linplot  <- qplot(medroi_fit$fitvals$Patlak_Ref, medroi_fit$fitvals$Patlak_ROI) + ggtitle('Medium') + xlab(patlak_xlab) + ylab(patlak_ylab)
  high_linplot <- qplot(highroi_fit$fitvals$Patlak_Ref, highroi_fit$fitvals$Patlak_ROI) + ggtitle('High') + xlab(patlak_xlab) + ylab(patlak_ylab)
  
  tstarInclFrames = 3:frames
  zeros = rep(0, length(tstarInclFrames))
  
  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  k_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  
  for(i in 1:length(tstarInclFrames)) {
    
    lowfit <- refPatlak(t_tac, reftac, lowroi, tstarIncludedFrames =  tstarInclFrames[i])
    medfit <- refPatlak(t_tac, reftac, medroi, tstarIncludedFrames =  tstarInclFrames[i])
    highfit <- refPatlak(t_tac, reftac, highroi, tstarIncludedFrames =  tstarInclFrames[i])
    
    r2_df$Low[i]    <- summary(lowfit$fit)$r.squared
    r2_df$Medium[i] <- summary(medfit$fit)$r.squared
    r2_df$High[i]   <- summary(highfit$fit)$r.squared
    
    maxperc_df$Low[i]    <- maxpercres(lowfit)
    maxperc_df$Medium[i] <- maxpercres(medfit)
    maxperc_df$High[i]   <- maxpercres(highfit)
    
    k_df$Low[i]    <- lowfit$par$K
    k_df$Medium[i] <- medfit$par$K
    k_df$High[i]   <- highfit$par$K
  }
  
  xlabel = 'Number of Included Frames' 
  ylab_r2 = expression(R^2)
  ylab_mp = 'Maximum Percentage Variance'
  
  
  # R Squared plots
  
  low_r2plot  <- ggplot(r2_df, aes(x=Frames, y = Low)) + geom_point() + scale_x_continuous(breaks=seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)
  med_r2plot  <- ggplot(r2_df, aes(x=Frames, y = Medium)) + geom_point() + scale_x_continuous(breaks=seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)
  high_r2plot <- ggplot(r2_df, aes(x=Frames, y = High)) + geom_point() + scale_x_continuous(breaks=seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0.99, 1)) + xlab(xlabel) + ylab(ylab_r2)
  
  # Max Percentage Variation Plots
  
  maxperc_df$inclmins <- rev( max(t_tac) - t_tac )[-c(1,2)]
  maxperc_df$tstar <- rev( t_tac)[-c(1,2)]
  
  low_mpplot  <- ggplot(maxperc_df, aes(x=Frames, y = Low)) + geom_point() + scale_x_continuous(breaks=seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text",x=3,y=20,label="t* Minutes",colour="red",size=3,hjust=0) + annotate("text",x=maxperc_df$Frames,y=maxperc_df$Low+1.4,label=round(maxperc_df$tstar,1),size=3,colour="red") + annotate("text",x=3,y=20-0.7,label="Included Minutes",colour="blue",size=3,hjust=0) + annotate("text",x=maxperc_df$Frames,y=maxperc_df$Low+0.7,label=round(maxperc_df$inclmins,1),size=3,colour="blue")
  med_mpplot  <- ggplot(maxperc_df, aes(x=Frames, y = Medium)) + geom_point() + scale_x_continuous(breaks=seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text",x=3,y=20,label="t* Minutes",colour="red",size=3,hjust=0) + annotate("text",x=maxperc_df$Frames,y=maxperc_df$Medium+1.4,label=round(maxperc_df$tstar,1),size=3,colour="red") + annotate("text",x=3,y=20-0.7,label="Included Minutes",colour="blue",size=3,hjust=0) + annotate("text",x=maxperc_df$Frames,y=maxperc_df$Medium+0.7,label=round(maxperc_df$inclmins,1),size=3,colour="blue")
  high_mpplot <- ggplot(maxperc_df, aes(x=Frames, y = High)) + geom_point() + scale_x_continuous(breaks=seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylim(c(0, 20)) + xlab(xlabel) + ylab(ylab_mp) + annotate("text",x=3,y=20,label="t* Minutes",colour="red",size=3,hjust=0) + annotate("text",x=maxperc_df$Frames,y=maxperc_df$High+1.4,label=round(maxperc_df$tstar,1),size=3,colour="red") + annotate("text",x=3,y=20-0.7,label="Included Minutes",colour="blue",size=3,hjust=0) + annotate("text",x=maxperc_df$Frames,y=maxperc_df$High+0.7,label=round(maxperc_df$inclmins,1),size=3,colour="blue")
  
  
  # TAC Plot
  
  tacplotdf <- data.frame(cbind( Time = lowroi_fit$tacs$Time , Reference = lowroi_fit$tacs$Reference, Low = lowroi_fit$tacs$Target ,  Medium = medroi_fit$tacs$Target , High = highroi_fit$tacs$Target ))
  tacplotdf <- tidyr::gather(tacplotdf, key = Region, value = Radioactivity, -Time)
  
  tacplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(tacplotdf$Region)))
  
  myColors <- RColorBrewer::brewer.pal(4,"Set1")
  names(myColors) <- levels(tacplotdf$Region)
  colScale <- scale_colour_manual(name = "Region",values = myColors)
  
  tacplot = ggplot(tacplotdf, aes(x=Time, y=Radioactivity, colour=Region)) + geom_point() + geom_line() + colScale
  
  
  # K Plot
  
  kplotdf <- tidyr::gather(k_df, key = Region, value = K, -Frames, -Time)
  kplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(kplotdf$Region)))
  
  kplot = ggplot(kplotdf, aes(x=Frames, y=K, colour=Region)) + geom_point() + geom_line() + scale_x_continuous(breaks=seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylab(expression(K[i])) + colScale
  
  
  # File Output
  
  jpeg(filename=paste0(filename,'_refPatlak.jpeg'),width=300,height=400,units="mm",res=600)
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(4, 3, heights = unit(c(1,1), "null"))))
  print(low_linplot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(med_linplot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(high_linplot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 3))
  print(low_r2plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(med_r2plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(high_r2plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 3))
  print(low_mpplot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  print(med_mpplot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 2))
  print(high_mpplot, vp = grid::viewport(layout.pos.row = 3, layout.pos.col = 3))
  print(tacplot, vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 1:2))
  print(kplot, vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 3))
  dev.off()
  
}