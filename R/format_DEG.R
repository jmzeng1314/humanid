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

  logFC_Cutof <- mean(abs(DEG$logFC)) + 2*sd(abs(DEG$logFC))
  table( abs(DEG$logFC) > logFC_Cutof & DEG$P.Value <0.05 )
  DEG$sigORnot <- ifelse(abs(DEG$logFC) > logFC_Cutof & DEG$P.Value <0.05 ,ifelse(DEG$logFC > logFC_Cutof,'UP','DOWN'),'NOT')
  table( DEG$sigORnot )

  png( paste0(studyID,'_volcanic.png') )
  plot(DEG$logFC,DEG$P.Value,main =logFC_Cutof )
  abline(h = 0.05,col='blue')
  abline(v=logFC_Cutof,col='red')
  abline(v=-logFC_Cutof,col='red')
  dev.off()

  file_allGeneList = paste0(studyID,'_allGeneList.txt')
  file_diffGeneList = paste0(studyID,'_diffGeneList.txt')
  file_upGeneList = paste0(studyID,'_upGeneList.txt')
  file_downGeneList = paste0(studyID,'_downGeneList.txt')

  write.table(unique(DEG$symbol),file_allGeneList,quote = F,row.names = F,col.names = F)
  write.table(unique(DEG[DEG$sigORnot != 'NOT','symbol']),file_diffGeneList,quote = F,row.names = F,col.names = F)
  write.table(unique(DEG[DEG$sigORnot == 'UP','symbol']), file_upGeneList ,quote = F,row.names = F,col.names = F)
  write.table(unique(DEG[DEG$sigORnot == 'DOWN','symbol']), file_downGeneList ,quote = F,row.names = F,col.names = F)

  batch_enrichment( file_diffGeneList, file_allGeneList ,studyID=paste0(studyID,'_diff') )
  batch_enrichment( file_upGeneList, file_allGeneList,studyID=paste0(studyID,'_UP') )
  batch_enrichment( file_downGeneList, file_allGeneList,studyID=paste0(studyID,'_DOWN') )


  logFC <- DEG$logFC;names(logFC) <- DEG$symbol

  keggEnrichTable=read.csv( paste0(studyID,'_diff_update_kegg.enrichment.csv') ,stringsAsFactors = F)
  keggEnrichTable=keggEnrichTable[keggEnrichTable$Pvalue<0.05,]
  diff_gene_list = read.table( file_diffGeneList ,stringsAsFactors = F)[,1]
  add_kegg_up_down_link(keggEnrichTable,logFC,diff_gene_list,studyID = paste0(studyID,'_diff'))

  keggEnrichTable=read.csv( paste0(studyID,'_UP_update_kegg.enrichment.csv') ,stringsAsFactors = F)
  keggEnrichTable=keggEnrichTable[keggEnrichTable$Pvalue<0.05,]
  diff_gene_list = read.table( file_diffGeneList ,stringsAsFactors = F)[,1]
  add_kegg_up_down_link(keggEnrichTable,logFC,diff_gene_list,studyID = paste0(studyID,'_UP'))

  keggEnrichTable=read.csv( paste0(studyID,'_DOWN_update_kegg.enrichment.csv') ,stringsAsFactors = F)
  keggEnrichTable=keggEnrichTable[keggEnrichTable$Pvalue<0.05,]
  diff_gene_list = read.table( file_diffGeneList ,stringsAsFactors = F)[,1]
  add_kegg_up_down_link(keggEnrichTable,logFC,diff_gene_list,studyID = paste0(studyID,'_DOWN'))



  tmpAnno <- geneAnno(unique(DEG$symbol));dim(tmpAnno)
  DEG <- merge(DEG,tmpAnno,by='symbol');dim(DEG)
  write.csv(DEG, paste0(studyID,'_DEG.csv'),row.names = F)


  rnk_file=paste0(studyID,'.rnk')
  sink( rnk_file )
  cat('#')
  sink()
  write.table(DEG[,c('symbol','logFC')] ,rnk_file,append = T,quote = F,row.names = F)


}
