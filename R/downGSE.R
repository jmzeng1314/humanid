#' download the data from GEO database
#'
#' Given a GSE study ID, just like GSE1009, this function will download the eSet object and write the expression matrix and phenotype information.
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1009
#' @param studyID A standard study ID in GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1009
#' @param destdir where to store the files just download.
#' @return A expresstionSet,S4 object in R, which contain: assayData,phenoData,featureData
#' @export
#' @keywords downGSE
#' @examples
#' #' downGSE('GSE1009')


downGSE <- function(studyID='GSE1009',destdir='.'){

  library(GEOquery)
  eSet <- getGEO(studyID, destdir=destdir,getGPL = F)

  exprSet=exprs(eSet[[1]])
  pdata=pData(eSet[[1]])

  write.csv(exprSet,paste0(studyID,'_exprSet.csv'))
  write.csv(pdata,paste0(studyID,'_metadata.csv'))
  return(eSet)

}
