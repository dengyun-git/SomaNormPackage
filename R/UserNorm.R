#' Takes in the raw RFUs data frame from SomaLogic adat file, and also user defined normalization workflow.
#' @param RawM the input RFUs data frame from SomaLogic Adat file.
#' @param Funlist user define a list: select from "HYBNORM,PLATESCALE,MIDNORM/MIDNORMcali/MIDNORMsamp,CALIBRATION.
#' @return MySoma the normalised RFUs after user defined normalization procedures.
#' @importFrom stats median
#' @export UserNorm

UserNorm <- function(Funlist,RawM){
  for (FunCounter in 1:length(Funlist)){
    f <- Funlist[[FunCounter]]
    MySoma = f(RawM)
    RawM = MySoma
  }
  return(MySoma)
}
