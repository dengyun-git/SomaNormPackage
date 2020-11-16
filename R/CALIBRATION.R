#' Calibration
#' Takes in RFUs. Output correponding RFUs after Calibration normalization.
#' @param RawM the input RFUs data frame from SomaLogic Adat file.
#' @return MySoma the output RFUs data frame after the Calibration normalization.
#' @importFrom stats median


CALIBRATION <- function(RawM){

  PlateIdUni = levels(factor(RawM$PlateId))

  Platelist = list()

  DatStartId <- which(colnames(RawM)=="CLI")+1  ###for calculation convenience, extract data zone only

  for (plateCounter in 1:length(PlateIdUni)){

    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])

    RawMS = RawM[PlateIdSg,] ### single plate

    idCaliborator = which(RawMS$SampleType=="Calibrator")

    datZone = RawMS[,DatStartId:ncol(RawMS)]

    CaliboratorM = datZone[idCaliborator,]

    CalibratorMedian = apply(CaliboratorM,2,median)

    CalSet=PlateScale_Reference/CalibratorMedian

    datZone2 = t(apply(datZone,1,function(x) x*CalSet))

    RawMSDone = cbind(RawMS[,1:DatStartId-1],datZone2)

    Platelist[[plateCounter]] = RawMSDone
  }

  MySoma = getMySoma(Platelist)

  return(MySoma)
}
