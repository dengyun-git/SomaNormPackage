#' Hybridizatioin normalization
#' Takes in RFUs. Output correponding RFUs after Hybridization normalization.
#' @param RawM the input RFUs data frame for Hybridizatioin normalisation
#' @return MySoma the output RFUs data frame after the hybridization normalization.
#' @importFrom stats median

HYBNORM <- function(RawM){

  PlateIdUni = levels(factor(RawM$PlateId))

  Platelist = list()

  DatStartId <- which(colnames(RawM)=="CLI")+1  ###for calculation convenience, extract data zone only

  HybId =  grep("HybControlElution",colnames(RawMS))

  for (plateCounter in 1:length(PlateIdUni)){

    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])

    RawMS = RawM[PlateIdSg,] ### single plate

    datZone = RawMS[,DatStartId:ncol(RawMS)] ###RFUs zone only

    rCmedian = matrix(apply(datZone[,HybId],2,median),nrow=1)

    HybNorm1 = t(apply(datZone[,HybId],1,function(x){rCmedian/x}))

    rRmedian = apply(HybNorm1,1,median)

    HybNorm = apply(datZone,2,function(x){x*rRmedian})

    RawMSDone = cbind(RawMS[,1:(DatStartId-1)],HybNorm)

    Platelist[[plateCounter]] = RawMSDone
  }

  MySoma = getMySoma(Platelist)

  return(MySoma)
}
