#' Median normalization
#' MIDNORM perform normalization across all the sample types included in the RawM$SampleType.
#' @param RawM the input RFUs data frame for median normalisation.
#' @return MySoma the output RFUs data frame after the Median normalisation.
#' @importFrom stats median

MIDNORM = function(RawM){ ###caliCase control which sample type to be applied MidNorm

  PlateIdUni = levels(factor(RawM$PlateId))

  Platelist = list()

  for (plateCounter in 1:length(PlateIdUni)){

    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])

    RawMS = RawM[PlateIdSg,] ### single plate

    idHyb = which(grepl("HybControlElution",colnames(RawMS))==TRUE)
    idNonHyb = which(!grepl("HybControlElution",colnames(RawMS))==TRUE)
    RawMS1 = RawMS[,idNonHyb] ###RawM1 to track "Ratio of Normalization Median to Sample Value"
    DatStartId = which(colnames(RawMS1) == "CLI") +1
    DatStartIdP = which(colnames(RawMS) == "CLI") +1

    sampType = levels(factor(RawMS$SampleType))

    for (sampTypeCounter in 1:length(sampType)){

      idSamp = which(RawMS1$SampleType == sampType[[sampTypeCounter]])

      datZoneP = RawMS1[,DatStartId:ncol(RawMS1)] ### datazone excludes "HybControlElution"

      datZone = RawMS1[idSamp,DatStartId:ncol(RawMS1)]  ###dataZone only include interested sample type

      if(length(idSamp)==1) {SampTypeRFU = matrix(datZone,nrow=1)
      SampTypeMedian = SampTypeRFU
      SampMedianNorm = as.matrix(SampTypeMedian/SampTypeRFU)}

      else {SampTypeRFU = datZone
      SampTypeMedian = apply(SampTypeRFU,2,median)
      SampMedianNorm = t(apply(SampTypeRFU,1,function(x){SampTypeMedian/x}))}

      Dilute = Dilution[which(Dilution!="0")] ###SampMedianNorm & Dilute: ncol(SampMedianNorm)=length(Dilute)
      uniqDilute = levels(factor(Dilute))

      for (idDilute in (1:length(uniqDilute))){

        DataDiluteID = which(Dilute==uniqDilute[idDilute])

        MedianDiluteSingle = matrix(apply(SampMedianNorm[,DataDiluteID],1,median),ncol=1)

        DataDiluteNormT = apply(datZone[,DataDiluteID],2,function(x){MedianDiluteSingle*x})

        if (idDilute==1){DataDiluteNorm = DataDiluteNormT}
        else{DataDiluteNorm = cbind(DataDiluteNorm,DataDiluteNormT)}
      }

      if (sampTypeCounter==1) {RawMStemp = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)}
      else{
        RawMStempT = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)
        RawMStemp = rbind(RawMStemp,RawMStempT)
      }
    }

    Platelist[[plateCounter]] = RawMStemp
  }

  MySomaTemp = getMySoma(Platelist)

  ###up to date MySomaTemp has different indexing from Raw. We make them consistent
  rowOrderName = rownames(RawM)
  rowOrder = vector(mode="numeric",length=nrow(MySomaTemp))
  for (j in 1:nrow(MySomaTemp)){
    rowOrder[j] = which(rownames(MySomaTemp) %in% rowOrderName[j])
  }

  colOrderName = colnames(RawM)
  colOrder = vector(mode="numeric",length=ncol(MySomaTemp))
  for (k in 1:ncol(MySomaTemp)){
    colOrder[k] = which(colnames(MySomaTemp) %in% colOrderName[k])
  }

  MySoma = MySomaTemp[rowOrder,colOrder]
  return(MySoma)
}


#' Median normalisation on calibrators
#' MIDNORMcali perform normalization on calibrators as shown in the RawM$SampleType.
#' @param RawM the input RFUs data frame for median normalisation.
#' @return MySoma the output RFUs data frame after the Median normalization on calibrators.
#'
MIDNORMcali = function(RawM){

  PlateIdUni = levels(factor(RawM$PlateId))

  Platelist = list()

  for (plateCounter in 1:length(PlateIdUni)){

    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])

    RawMS = RawM[PlateIdSg,] ### single plate

    idHyb = which(grepl("HybControlElution",colnames(RawMS))==TRUE)
    idNonHyb = which(!grepl("HybControlElution",colnames(RawMS))==TRUE)
    RawMS1 = RawMS[,idNonHyb] ###RawM1 to track "Ratio of Normalization Median to Sample Value"
    DatStartId = which(colnames(RawMS1) == "CLI") +1
    DatStartIdP = which(colnames(RawMS) == "CLI") +1

    sampType = "Calibrator"

    for (sampTypeCounter in 1:length(sampType)){

      idSamp = which(RawMS1$SampleType == sampType[[sampTypeCounter]])

      datZoneP = RawMS1[,DatStartId:ncol(RawMS1)] ### datazone excludes "HybControlElution"

      datZone = RawMS1[idSamp,DatStartId:ncol(RawMS1)]  ###dataZone only include interested sample type

      if(length(idSamp)==1) {SampTypeRFU = matrix(datZone,nrow=1)
      SampTypeMedian = SampTypeRFU
      SampMedianNorm = as.matrix(SampTypeMedian/SampTypeRFU)}

      else {SampTypeRFU = datZone
      SampTypeMedian = apply(SampTypeRFU,2,median)
      SampMedianNorm = t(apply(SampTypeRFU,1,function(x){SampTypeMedian/x}))}

      Dilute = Dilution[which(Dilution!="0")] ###SampMedianNorm & Dilute: ncol(SampMedianNorm)=length(Dilute)
      uniqDilute = levels(factor(Dilute))

      for (idDilute in (1:length(uniqDilute))){

        DataDiluteID = which(Dilute==uniqDilute[idDilute])

        MedianDiluteSingle = matrix(apply(SampMedianNorm[,DataDiluteID],1,median),ncol=1)

        DataDiluteNormT = apply(datZone[,DataDiluteID],2,function(x){MedianDiluteSingle*x})

        if (idDilute==1){DataDiluteNorm = DataDiluteNormT}
        else{DataDiluteNorm = cbind(DataDiluteNorm,DataDiluteNormT)}
      }

      if (sampTypeCounter==1) {RawMStemp = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)}
      else{
        RawMStempT = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)
        RawMStemp = rbind(RawMStemp,RawMStempT)
      }
      RawMStemp2 = rbind(RawMStemp,RawMS[which(RawMS$SampleType!="Calibrator"),])
    }

    Platelist[[plateCounter]] = RawMStemp2
  }

  MySomaTemp = getMySoma(Platelist)

  ###up to date MySomaTemp has different indexing from Raw. We make them consistent
  rowOrderName = rownames(RawM)
  rowOrder = vector(mode="numeric",length=nrow(MySomaTemp))
  for (j in 1:nrow(MySomaTemp)){
    rowOrder[j] = which(rownames(MySomaTemp) %in% rowOrderName[j])
  }

  colOrderName = colnames(RawM)
  colOrder = vector(mode="numeric",length=ncol(MySomaTemp))
  for (k in 1:ncol(MySomaTemp)){
    colOrder[k] = which(colnames(MySomaTemp) %in% colOrderName[k])
  }

  MySoma = MySomaTemp[rowOrder,colOrder]
  return(MySoma)
}


#' Median normalisation on samples excluding Calibrators
#' @param RawM the input RFUs data frame for median normalization
#' @return MySoma the output RFUs data frame after the Median normalization on samples
MIDNORMsamp = function(RawM){

  PlateIdUni = levels(factor(RawM$PlateId))

  Platelist = list()

  for (plateCounter in 1:length(PlateIdUni)){

    PlateIdSg = which(RawM$PlateId == PlateIdUni[plateCounter])

    RawMS = RawM[PlateIdSg,] ### single plate

    idHyb = which(grepl("HybControlElution",colnames(RawMS))==TRUE)
    idNonHyb = which(!grepl("HybControlElution",colnames(RawMS))==TRUE)
    RawMS1 = RawMS[,idNonHyb] ###RawM1 to track "Ratio of Normalization Median to Sample Value"
    DatStartId = which(colnames(RawMS1) == "CLI") +1
    DatStartIdP = which(colnames(RawMS) == "CLI") +1

    sampTypePre = levels(factor(RawMS$SampleType))
    sampType = which(sampTypePre!="Calibrator")

    for (sampTypeCounter in 1:length(sampType)){

      idSamp =grep(sampType[sampTypeCounter],RawMS1$SampleType)

      datZoneP = RawMS1[,DatStartId:ncol(RawMS1)] ### datazone excludes "HybControlElution"

      datZone = RawMS1[idSamp,DatStartId:ncol(RawMS1)]  ###dataZone only include interested sample type

      if(length(idSamp)==1) {SampTypeRFU = matrix(as.matrix(datZone),nrow=1)
      SampTypeMedian = SampTypeRFU
      SampMedianNorm = as.matrix(SampTypeMedian/SampTypeRFU)}

      else {SampTypeRFU = datZone
      SampTypeMedian = apply(SampTypeRFU,2,median)
      SampMedianNorm = t(apply(SampTypeRFU,1,function(x){SampTypeMedian/x}))}

      Dilute = Dilution[which(Dilution!="0")] ###SampMedianNorm & Dilute: ncol(SampMedianNorm)=length(Dilute)
      uniqDilute = levels(factor(Dilute))

      for (idDilute in (1:length(uniqDilute))){

        DataDiluteID = which(Dilute==uniqDilute[idDilute])

        MedianDiluteSingle = matrix(apply(SampMedianNorm[,DataDiluteID],1,median),ncol=1)

        DataDiluteNormT = apply(datZone[,DataDiluteID],2,function(x){MedianDiluteSingle*x})

        if (idDilute==1){DataDiluteNorm = DataDiluteNormT}
        else{DataDiluteNorm = cbind(DataDiluteNorm,DataDiluteNormT)}
      }

      if (sampTypeCounter==1) {RawMStemp = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)}
      else{
        RawMStempT = cbind(RawMS[idSamp,1:(DatStartIdP-1)],RawMS[idSamp,idHyb],DataDiluteNorm)
        RawMStemp = rbind(RawMStemp,RawMStempT)
      }
      RawMStemp2 = rbind(RawMStemp,RawMS[which(RawMS$SampleType =="Calibrator"),])
    }

    Platelist[[plateCounter]] = RawMStemp2
  }

  MySomaTemp = getMySoma(Platelist)

  ###up to date MySomaTemp has different indexing from Raw. We make them consistent
  rowOrderName = rownames(RawM)
  rowOrder = vector(mode="numeric",length=nrow(MySomaTemp))
  for (j in 1:nrow(MySomaTemp)){
    rowOrder[j] = which(rownames(MySomaTemp) %in% rowOrderName[j])
  }

  colOrderName = colnames(RawM)
  colOrder = vector(mode="numeric",length=ncol(MySomaTemp))
  for (k in 1:ncol(MySomaTemp)){
    colOrder[k] = which(colnames(MySomaTemp) %in% colOrderName[k])
  }

  MySoma = MySomaTemp[rowOrder,colOrder]
  return(MySoma)
}
