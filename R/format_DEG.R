#' write some files from the DEG results from topTable(limma)
#'
#'
#' @param DEG  DEG=topTable(fit,coef=2,adjust='BH')
#' @param studyID A standard study ID
#' @return write some files
#' @export
#' @keywords format_DEG
#' @examples
#' #'  format_DEG(DEG)
#'
#'
format_DEG <- function(DEG,studyID='test'){
  DEG$symbol=rownames(DEG)
  logFC_Cutof <- mean(abs(DEG$logFC)) + 2*sd(abs(DEG$logFC))
  table( abs(DEG$logFC) > logFC_Cutof & DEG$P.Value <0.05 )
  DEG$sigORnot <- abs(DEG$logFC) > logFC_Cutof & DEG$P.Value <0.05

  write.table(unique(DEG$symbol),paste0(studyID,'_allGeneList.txt'),quote = F,row.names = F,col.names = F)
  write.table(unique(DEG[DEG$sigORnot,'symbol']),paste0(studyID,'_diffGeneList.txt'),quote = F,row.names = F,col.names = F)

  rnk_file=paste0(studyID,'.rnk')
  sink( rnk_file )
  cat('#')
  sink()
  write.table(DEG[,c('symbol','logFC')] ,rnk_file,append = T,quote = F,row.names = F)


  tmpAnno <- geneAnno(DEG$symbol)
  DEG <- merge(DEG,tmpAnno,by='symbol')
  write.csv(DEG, paste0(studyID,'_DEG.csv'),row.names = F)

}
