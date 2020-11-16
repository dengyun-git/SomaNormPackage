#' Initialise required metadata for normalisation methods
#' @param inputfile1 SomaLogic Adat file, which contains the raw RFUs
#' @param inputfile2 SomaLogic Adat file, which contains metadata.Only external reference is required from this file.
#' @return RawM the raw RFUs data frame for normalisation.
#' @importFrom  SomaDataIO read.adat
#' @export initQCnorm

initQCnorm = function(inputfile1,inputfile2){

  RawM <- SomaDataIO::read.adat(inputfile1)

  ### read in Col^MetaTable
  ###use inputfile1, there is no PlateScale_Reference information, so we use inputfile2 considering all the row and column indices are consistent

  con1 = file(inputfile2, "r")

  while(TRUE) {
    sline = readLines(con1, n=1)

    if(length(sline) == 0){
      print("Meta data not intact")
      break}

    slineL = strsplit(sline,'\t')[[1]]

    if(length(slineL) > 29){

      if(slineL[30]=="Dilution"){Dilution <<- slineL[31:length(slineL)]}
      if(slineL[30]=="PlateScale_Reference"){
        PlateScale_Reference <<- as.numeric(slineL[31:length(slineL)])
        break}
      ### 30 is the beginning field of the ^COL_DATA;column names of meta data = SeqId +Target
    }
  }
  close(con1)

  return(RawM)
}
