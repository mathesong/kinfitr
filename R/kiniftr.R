#' @importFrom stats AIC as.formula coef fitted lag lm median na.omit pgamma predict residuals start time weights
#' @importFrom utils data head modifyList read.delim tail
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics lines
#' @importFrom rlang :=
NULL

# Suppress R CMD check notes about NSE variables used in dplyr/ggplot2 functions
utils::globalVariables(c(
  ".", ".x", "AIF", "AIFmodel", "BP", "BPR", "Blood", "Data", "Equilibrium", 
  "Exclude", "Fitted", "Frames", "High", "K", "Label", "Logan_Plasma", 
  "Logan_ROI", "Logan_ref", "Logan_roi", "Low", "Measurement", "Medium", 
  "Method", "Outcome", "ParentFraction", "Patlak_Plasma", "Patlak_ROI", 
  "Patlak_ref", "Patlak_roi", "Plasma", "Pred", "RSS", "RSSw", "Radioactivity", 
  "Region", "TAC", "Term1_DV", "Time", "Unit_AIF", "Unit_time", "Value", 
  "Vnd", "Vt", "Weights", "acq", "activity", "activity_1ex", "activity_2ex", 
  "all_of", "attribute", "data", "dotsize", "dur", "everything", "extension", 
  "filedata", "fit", "inpshift", "jsondata", "keep", "lag", "lead", "log_RSS", 
  "measplasma", "measurement", "metabolite_parent_fraction", "parentFraction", 
  "path", "path_absolute", "path_relative", "plasma_radioactivity", 
  "plasma_uncor", "prop_success", "rec", "recording", "run", "sampleDuration", 
  "sampleStartTime", "ses", "start", "starts_with", "success", "tactimes", 
  "task", "time", "trc", "tsvdata", "update_blooddata_bids", "value", 
  "whole_blood_radioactivity", "ROI.fitted", "ROI.measured", "ROI_measured", 
  "Reference", "Target.fitted", "Target.measured", "Measure"
))
