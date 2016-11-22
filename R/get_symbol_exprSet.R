#' extract the expression matrix from a eSet object
#'
#' not only extract the expression matrix, but change the probeset ID to gene HUGO symbol ,and remove the duplicated genes .
#'
#' @param eSet  A expresstionSet,S4 object in R, which contain: assayData,phenoData,featureData
#' @param platformDB which microarray platform does the eSet use ?
#' @return expression matrix, rowname is the gene HUGO symbol
#' @export
#' @keywords get_symbol_exprSet
#' @examples
#' #' exprSet <- get_symbol_exprSet(K27M_WT_eSet)
#'
#'
get_symbol_exprSet <- function(eSet,platformDB='hgu133plus2.db'){
  ## you must make sure that the eSet is a standard format, eg: read the cel files by Affy package
  library(platformDB, character.only=TRUE)
  library(annotate)
  probeset <- featureNames( eSet )
  exprSet <- exprs( eSet)
  #EGID <- as.numeric(lookUp(probeset, platformDB, "ENTREZID"))
  SYMBOL <-  lookUp(probeset, platformDB, "SYMBOL")
  length(unlist(SYMBOL));dim(exprSet)

  exprSet=as.data.frame(exprSet)
  a=cbind(as.character(unlist(SYMBOL)),exprSet)

  exprSet=rmDupID(a)

  keepProbe <- apply(exprSet, 1, function(x) all( x >1 ) )
  exprSet <- exprSet[keepProbe,];dim(exprSet)

  if( mean(rowMeans( exprSet ,na.rm = T),na.rm = T) >20)
          exprSet=log2(exprSet) ## based on 2

  #boxplot(exprSet,las=2)
  #exprSet <- na.omit(exprSet)
  return(exprSet)

}
