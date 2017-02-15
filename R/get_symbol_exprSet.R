#' extract the expression matrix from a eSet object
#'
#' not only extract the expression matrix, but change the probeset ID to gene HUGO symbol ,and remove the duplicated genes .
#'
#' @param eSet  A expresstionSet,S4 object in R, which contain: assayData,phenoData,featureData
#' @param platformDB which microarray platform does the eSet use,such as  illuminaHumanv4.db, hgu133plus2.db
#' @param filter whether to filter genes show very low expression values, default:TRUE
#' @return expression matrix, rowname is the gene HUGO symbol
#' @export
#' @keywords get_symbol_exprSet
#' @examples
#' #' exprSet <- get_symbol_exprSet(K27M_WT_eSet)
#'
#'
get_symbol_exprSet <- function(eSet,platformDB='hgu133plus2.db',filter=TRUE){
  ## you must make sure that the eSet is a standard format, eg: read the cel files by Affy package
  library(platformDB, character.only=TRUE)
  library(annotate)
  probeset <- featureNames( eSet )
  exprSet <- exprs( eSet)
  exprSet <- na.omit(exprSet )

  probe2symbol_df <- toTable(get(paste0(sub('.db','',platformDB),'SYMBOL') ))

  exprSet=as.data.frame(exprSet)
  exprSet$probe_id = rownames(exprSet)

  tmp <- merge(probe2symbol_df,exprSet,by='probe_id');dim(tmp)
  tmp <- tmp[,-1]
  exprSet_rmdup=rmDupID(tmp);dim(exprSet_rmdup)
  if(filter){
    keepProbe <- apply(exprSet_rmdup, 1, function(x) all( x >1 ) )
    exprSet_rmdup <- exprSet_rmdup[as.logical(keepProbe),];dim(exprSet_rmdup)

  }

  if( mean(rowMeans( exprSet_rmdup ,na.rm = T),na.rm = T) >20)
    exprSet_rmdup=log2(exprSet_rmdup) ## based on 2

  #boxplot(exprSet_rmdup,las=2)
  #exprSet_rmdup <- na.omit(exprSet_rmdup)
  return(exprSet_rmdup)

}
