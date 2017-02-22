#' Ichise Multilinear Analysis 1
#'
#' Function to fit the MA1 of Ichise et al (2002) to data.
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a 
#' zero. If a time zero frame is not included, it will be added.
#' @param tac Numeric vector of radioactivity concentrations in the target tissue for each frame. We include zero at time 
#' zero: if not included, it is added.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param tstarIncludedFrames The number of frames to be used in the regression model, i.e. the number of frames for which 
#' the function is linear after pseudo-equilibrium is reached. This is a count from the end of the measurement, so a value of 
#' 10 means that last 10 frames will be used. This value can be estimated using \code{ma1_tstar}.
#' @param weights Optional. Numeric vector of the weights assigned to each frame in the fitting. We include zero at time zero: 
#' if not included, it is added. If not specified, uniform weights will be used.
#' @param inpshift Optional. The number of minutes by which to shift the timing of the input data frame forwards or backwards.
#' If not specified, this will be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param vB Optional. The blood volume fraction.  If not specified, this will be ignored and assumed to be 0%. If specified, it
#' will be corrected for prior to parameter estimation using the following equation: 
#' \deqn{C_{T}(t) = \frac{C_{Measured}(t) - vB\times C_{B}(t)}{1-vB}}
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20). 
#' This is to assess time stability.
#'
#' @return A list with a data frame of the fitted parameters \code{out$par}, the model fit object \code{out$fit}, 
#' a dataframe containing the TACs of the data \code{out$tacs}, a dataframe containing the fitted values \code{out$fitvals},
#' the blood input data frame after time shifting \code{input}, a vector of the weights \code{out$weights}, 
#' the inpshift value used \code{inpshift}, the specified vB value \code{out$vB} and the specified tstarIncludedFrames 
#' value \code{out$tstarIncludedFrames}.
#'
#' @examples
#' ma1(t_tac, tac, input, 10, weights, inpshift = onetcmout$par$inpshift)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#' 
#' @references Ichise M, Toyama H, Innis RB, Carson RE. Strategies to improve neuroreceptor parameter estimation by linear regression analysis. Journal of Cerebral Blood Flow & Metabolism. 2002 Oct 1;22(10):1271-81.
#'
#' @export

ma1 <- function(t_tac, tac, input, tstarIncludedFrames, weights, inpshift = 0, vB = 0, frameStartEnd) {
  
  
  # Tidying
  
  tidyinput <- tidyinput_art(t_tac, tac, weights, frameStartEnd)
  
  t_tac   <- tidyinput$t_tac
  tac  <- tidyinput$tac
  weights <- tidyinput$weights

  newvals = shift_timings(t_tac = t_tac,
                          tac = tac,
                          input = input, 
                          inpshift = inpshift )
  
  t_tac <- newvals$t_tac
  tac <- newvals$tac
  
  t_inp <- newvals$input$Time
  blood <- newvals$input$blood
  plasma <- newvals$input$plasma
  parentfrac <- newvals$input$parentfrac
  
  # Parameters
  
  corrplasma <- newvals$input$plasma * newvals$input$parentfrac
  
  interptime <- newvals$input$Time
  
  i_tac <- pracma::interp1(t_tac, tac, interptime, method="linear")
  
    # Blood Volume Correction (nothing happens if vB = 0)
    i_tac <- ( i_tac - vB*blood ) / ( 1 - vB )
  
  term1 = as.numeric( pracma::cumtrapz(interptime, corrplasma) )
  term2 = as.numeric( pracma::cumtrapz(interptime, i_tac) )
  
  term1 <- pracma::interp1(interptime, term1, t_tac, method="linear")
  term2 <- pracma::interp1(interptime, term2, t_tac, method="linear")
  
  tac_equil = tail(tac, tstarIncludedFrames)
  t_tac_equil = tail(t_tac, tstarIncludedFrames)
  term1_equil = tail(term1, tstarIncludedFrames)
  term2_equil = tail(term2, tstarIncludedFrames)
  weights_equil   = tail(weights, tstarIncludedFrames)
  
  
  # Solution
  
  ma1_model <- lm(tac_equil ~ term1_equil + term2_equil - 1, weights = weights_equil )
  
  # Output
  
  Vt = as.numeric( -ma1_model$coefficients[1] / ma1_model$coefficients[2])
  
  par = as.data.frame(list(Vt = Vt))
  fit = ma1_model
  
  tacs = data.frame(Time = t_tac, Target = tac)
  
  fitvals = data.frame(Time = t_tac_equil, Target = tac_equil, Term1 = term1_equil, Term2 = term2_equil, 
                       Target_fitted = as.numeric(predict(ma1_model)), Weights = weights_equil )
  
  input <- newvals$input
  
  out <- list(par = par, fit = ma1_model, tacs = tacs, fitvals = fitvals, 
              input = input, weights = weights, inpshift = inpshift, vB = vB, 
              tstarIncludedFrames = tstarIncludedFrames, model='ma1')
  return(out)
  
}

#' Plot: Ichise Multilinear Analysis 1
#'
#' Function to visualise the fit of the MA1 model to data.
#'
#' @param ma1out The output object of the MA1 fitting procedure.
#' @param roiname Optional. The name of the Target Region to see it on the plot.
#'
#' @return A ggplot2 object of the plot.
#'
#' @examples
#' plot_ma1fit(ma1out)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#' 
#' @import ggplot2
#'
#' @export

plot_ma1fit <- function(ma1out, roiname = NULL) {
 
  measured <- data.frame(Time = ma1out$tacs$Time,
                         Target.measured = ma1out$tacs$Target,
                         Weights = ma1out$weights)
  
  fitted <- data.frame(Time = ma1out$fitvals$Time,
                       Target.fitted = ma1out$fitvals$Target_fitted,
                       Weights = ma1out$fitvals$Weights )
  
  if(is.null(roiname)) {roiname = 'ROI'}
  
  measured = plyr::rename(measured, c('Target.measured' = paste0(roiname, '.measured')) )
                          
  fitted = plyr::rename(fitted, c('Target.fitted' = paste0(roiname, '.fitted')) )
  
  tidymeasured <- tidyr::gather(measured, key=Region, value=Radioactivity,
                                -Time, -Weights, factor_key = F)
  
  tidyfitted <- tidyr::gather(fitted, key=Region, value=Radioactivity,
                              -Time, -Weights, factor_key = F)
  
  Region <- forcats::fct_inorder(factor(c(tidymeasured$Region, tidyfitted$Region)) )
  
  myColors <- RColorBrewer::brewer.pal(3,"Set1")
  names(myColors) <- levels(Region)
  colScale <- scale_colour_manual(name = "Region",values = myColors)
  
  outplot <- ggplot(tidymeasured, aes(x=Time, y=Radioactivity, colour=Region)) + 
    geom_point(data=tidymeasured, aes(shape='a', size=Weights)) + 
    geom_line(data=tidyfitted) + colScale +
    guides(shape=FALSE, color=guide_legend(order=1)) + scale_size(range=c(1,3))
  
  return(outplot)
   
}


#' Tstar Finder: Ichise Multilinear Analysis 1
#'
#' Function to identify where t* is for MA1.
#'
#'
#' @param t_tac Numeric vector of times for each frame in minutes. We use the time halfway through the frame as well as a 
#' zero. If a time zero frame is not included, it will be added.
#' @param lowroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with low binding.
#' @param medroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with medium binding.
#' @param highroi Numeric vector of radioactivity concentrations in a target tissue for each frame. This should be from a ROI with high binding.
#' @param input Data frame containing the blood, plasma, and parent fraction concentrations over time.  This can be generated
#' using the \code{blood_interp} function.
#' @param filename The name of the output image: filename_ma1.jpeg
#' @param inpshift Optional. The number of minutes by which to shift the timing of the input data frame forwards or backwards.
#' If not specified, this will be set to 0. This can be fitted using 1TCM or 2TCM.
#' @param frameStartEnd Optional: This allows one to specify the beginning and final frame to use for modelling, e.g. c(1,20). 
#' This is to assess time stability.
#' @param gridbreaks Optional. The size of the grid in the plots. Default: 2.
#'
#' @return Saves a jpeg of the plots as filename_ma1.jpeg
#'
#' @examples
#' ma1_tstar(t_tac, lowroi, medroi, highroi, input, filename='demonstration', inpshift = onetcmout$par$inpshift, frameStartEnd, gridbreaks=4)
#'
#' @author Granville J Matheson, \email{mathesong@@gmail.com}
#' 
#' @import ggplot2
#'
#' @export


ma1_tstar <- function(t_tac, lowroi, medroi, highroi, input, filename, inpshift = 0, frameStartEnd, gridbreaks=2) {
  
  frames = length(t_tac)
  lowroi_fit  <- ma1(t_tac, lowroi, input, tstarIncludedFrames = frames, inpshift = inpshift, frameStartEnd = frameStartEnd)
  medroi_fit  <- ma1(t_tac, medroi, input, tstarIncludedFrames = frames, inpshift = inpshift, frameStartEnd = frameStartEnd)
  highroi_fit <- ma1(t_tac, highroi, input, tstarIncludedFrames = frames, inpshift = inpshift, frameStartEnd = frameStartEnd)
  
  low_linplot  <- plot_ma1fit(lowroi_fit) + ggtitle('Low') + ylim(0, max(lowroi_fit$tacs$Target*1.1))  + theme(legend.position="none")
  med_linplot  <- plot_ma1fit(medroi_fit) + ggtitle('Medium') + ylim(0, max(medroi_fit$tacs$Target*1.1))  + theme(legend.position="none")
  high_linplot <- plot_ma1fit(highroi_fit) + ggtitle('High') + ylim(0, max(highroi_fit$tacs$Target*1.1))  + theme(legend.position="none")
  
  tstarInclFrames = 3:frames
  zeros = rep(0, length(tstarInclFrames))
  
  r2_df <- data.frame(Frames = tstarInclFrames, Low = zeros, Medium = zeros, High = zeros)
  maxperc_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  vt_df <- data.frame(Frames = tstarInclFrames, Time = t_tac[ tstarInclFrames ], Low = zeros, Medium = zeros, High = zeros)
  
  for(i in 1:length(tstarInclFrames)) {
    
    lowfit  <- ma1(t_tac, lowroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, frameStartEnd = frameStartEnd)
    medfit  <- ma1(t_tac, medroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, frameStartEnd = frameStartEnd)
    highfit <- ma1(t_tac, highroi, input, tstarIncludedFrames = tstarInclFrames[i], inpshift = inpshift, frameStartEnd = frameStartEnd)
    
    r2_df$Low[i]    <- summary(lowfit$fit)$r.squared
    r2_df$Medium[i] <- summary(medfit$fit)$r.squared
    r2_df$High[i]   <- summary(highfit$fit)$r.squared
    
    maxperc_df$Low[i]    <- maxpercres(lowfit)
    maxperc_df$Medium[i] <- maxpercres(medfit)
    maxperc_df$High[i]   <- maxpercres(highfit)
    
    vt_df$Low[i]    <- lowfit$par$Vt
    vt_df$Medium[i] <- medfit$par$Vt
    vt_df$High[i]   <- highfit$par$Vt
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
  
  tacplotdf <- data.frame(cbind( Time = lowroi_fit$tacs$Time , Low = lowroi_fit$tacs$Target ,  Medium = medroi_fit$tacs$Target , High = highroi_fit$tacs$Target ))
  tacplotdf <- tidyr::gather(tacplotdf, key = Region, value = Radioactivity, -Time)
  
  tacplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(tacplotdf$Region)))
  
  myColors <- RColorBrewer::brewer.pal(4,"Set1")
  names(myColors) <- levels(tacplotdf$Region)
  colScale <- scale_colour_manual(name = "Region",values = myColors)
  
  tacplot = ggplot(tacplotdf, aes(x=Time, y=Radioactivity, colour=Region)) + geom_point() + geom_line() + colScale
  
  
  # Vt Plot
  
  vtplotdf <- tidyr::gather(vt_df, key = Region, value = Vt, -Frames, -Time)
  vtplotdf$Region <- forcats::fct_rev(forcats::fct_inorder(factor(vtplotdf$Region)))
  
  ylimits = c( min(vtplotdf$Vt) , max(vtplotdf$Vt) )
  if( max(vtplotdf$Vt) > 20 || min(vtplotdf$Vt) < 0 ) { ylimits = c( 0 , 20 ) }
  
  vtplot = ggplot(vtplotdf, aes(x=Frames, y=Vt, colour=Region)) + geom_point() + geom_line() + scale_x_continuous(breaks=seq(min(tstarInclFrames), max(tstarInclFrames), by = gridbreaks)) + ylab(expression(V[T])) + ylim(ylimits)  + colScale
  
  
  # File Output
  
  jpeg(filename=paste0(filename,'_ma1.jpeg'),width=300,height=400,units="mm",res=600)
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
  print(vtplot, vp = grid::viewport(layout.pos.row = 4, layout.pos.col = 3))
  dev.off()
  
}