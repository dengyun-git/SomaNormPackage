#' Combine platewise normalised data frame into one large data frame.
#' Takes in a list with elements of normalized data frame per plate. Output a whole data frame combining the element of Platelist.
#' @param Platelist the input RFUs data frame per plate.
#' @return MySoma the whole RFUs data frame containing all the plates


getMySoma <- function(Platelist){
  MySoma = data.frame(Platelist[[1]]) ###MySoma: the whole dataframe for all the plates
  for (plateC in 2:length(Platelist)){
    MySomaT = data.frame(Platelist[[plateC]])
    MySoma = rbind(MySoma,MySomaT)
  }
  return(MySoma)
}
