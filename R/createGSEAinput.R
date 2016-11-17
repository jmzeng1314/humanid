#' create gct file and cls file
#'
#' create gct file and cls file for GSEA according the expression matrix and group information
#'
#' @param studyID A standard study ID in GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49822
#' @param exprSet  a expression matrix, which columns are sample,rows are HUGO gene symbols.
#' @param group_list a vector,as long as the col number for the expression matrix,which describe the group for the samples in exprSet
#' @param destdir where to store the files just download.
#' @return A expresstionSet,S4 object in R, which contain: assayData,phenoData,featureData
#' @export
#' @keywords downGSE
#' @examples
#' #' createGSEAinput('GSE1009')


createGSEAinput <- function(studyID='GSE1009',exprSet=example_exprSet ,group_list ,destdir='.'){
  # sink("outfile.txt")
  # cat("hello")
  # cat("\n")
  # cat("world")
  # sink()
  gct_file=paste0(studyID,'.gct')
  sink( gct_file )
  cat("#1.2\n")
  cat(paste0( nrow(exprSet) ,"\t",length(unique(group_list)) ,"\n") )
  sink()
  gct_out <- cbind(symbol=rownames(exprSet),description='na',exprSet)
  write.table(gct_out,gct_file,append = T,quote = F,row.names = F)


  cls_file=paste0(studyID,'.cls')
  sink( cls_file )
  cat(paste0( length(group_list)  ,"\t",length(unique(group_list)) ,"\t1\n") )
  cat(paste0("#",paste(unique(group_list),collapse = '\t') ,"\n"))
  cat(paste(group_list,collapse = '\t'))
  sink()
}



