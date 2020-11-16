#' Plate scaling
#' @param RawM the input RFUs data frame for plate scaling.
#' @return MySoma the output RFUs data frame after the plate scaling normalization.
#' @importFrom stats median

PLATESCALE <- function(RawM){

  PlateIdUni = levels(factor(RawM$PlateId))

  Platelist = list()

  DatStartId <- which(colnames(RawM)=="CLI")+1  ###for calculation convenience, extract data zone only

  for (plateCounter in 1:length(PlateIdUni)){

    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])

    RawMS = RawM[PlateIdSg,] ### single plate

    idCaliborator = which(RawMS$SampleType=="Calibrator")

    datZone = RawMS[,DatStartId:ncol(RawMS)]

    CaliboratorM = datZone[idCaliborator,]

    calibratorMedian = apply(CaliboratorM,2,median)

    PlateScaleRatio = PlateScale_Reference/calibratorMedian

    PlateScaleScalar = median(PlateScaleRatio)

    datZone2 = datZone*PlateScaleScalar

    RawMSDone = cbind(RawMS[,1:(DatStartId-1)],datZone2)

    Platelist[[plateCounter]] = RawMSDone
  }

  MySoma = getMySoma(Platelist)

  return(MySoma)
}

