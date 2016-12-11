#' enrichment for GO/KEGG
#'
#'
#' @param diff_gene_file   just one column, gene HUGO symbol or entrez ID , no need for column name.
#' @param all_genes_file   the same to the diff_gene_file, contain more genes than it.
#' @param studyID A standard study ID
#' @param GOstats whether use GOstats or not , defaults:F
#' @return write enrichment results into files.
#' @export
#' @keywords batch_enrichment
#' @examples
#' #'  batch_enrichment('diff_gene_file.txt','all_genes_file.txt')
#'
#'
batch_enrichment <- function(diff_gene_file,all_genes_file,studyID='test',GOstats=F ){

  diff_gene_list = read.table(diff_gene_file,stringsAsFactors = F)[,1]
  all_genes_list = read.table(all_genes_file,stringsAsFactors = F)[,1]
  diff_gene_list =  unique(geneAnno(diff_gene_list)$gene_id)
  all_genes_list =  unique(geneAnno(all_genes_list)$gene_id)


  kegg_result=hyperGtest_jimmy(GeneID2kegg,kegg2GeneID,diff_gene_list,all_genes_list)
  kegg_result=as.data.frame(kegg_result,stringsAsFactors = F)
  kegg_result$pathway_name=kegg2name[match(kegg_result[,1],kegg2name[,'pathway_id']),'pathway_name']
  kegg_result=kegg_result[order(as.numeric(kegg_result[,2])),]
  kegg_result$PathwayID = paste0('hsa:', kegg_result$PathwayID )
  kegg_result$p.adjust = p.adjust(kegg_result$Pvalue,method = 'BH')
  write.csv(kegg_result,paste0(studyID,"_update_kegg.enrichment.csv"),row.names = F)

  if(GOstats){
    annotationPKG='org.Hs.eg.db'
    suppressMessages(library(GO.db))
    suppressMessages(library(KEGG.db))
    suppressMessages(library(GOstats))
    suppressMessages(library(org.Hs.eg.db))

    if(T){
      GOES = c('BP','CC', 'MF');
      for (ontology in GOES) {
        GO.hyperG.params = new("GOHyperGParams", geneIds=diff_gene_list, universeGeneIds=all_genes_list, annotation=annotationPKG,
                               ontology=ontology, pvalueCutoff=1, conditional = FALSE, testDirection = "over")
        GO.hyperG.results = hyperGTest(GO.hyperG.params);
        outHTMLname=paste("GO_",ontology,".enrichment.html",sep="")
        #htmlReport(GO.hyperG.results, file=outHTMLname, summary.args=list("htmlLinks"=TRUE))
        GO.hyperG.matrix=summary(GO.hyperG.results)
        outMatrixName=paste("GO.",ontology,".hyperG.summary.csv",sep="")
        write.csv(GO.hyperG.matrix,outMatrixName )

        outHTMLname=paste("GO.",ontology,".hyperG.summary.html",sep="")
        htmlReport(GO.hyperG.results, file=outHTMLname, summary.args=list("htmlLinks"=TRUE))

      }
      #options(digits=4);
      hyperG.params = new("KEGGHyperGParams", geneIds=diff_gene_list, universeGeneIds=all_genes_list, annotation=annotationPKG,
                          categoryName="KEGG", pvalueCutoff=1, testDirection = "over")
      KEGG.hyperG.results = hyperGTest(hyperG.params);
      outHTMLname="kegg.enrichment.html"
      htmlReport(KEGG.hyperG.results, file=outHTMLname, summary.args=list("htmlLinks"=TRUE))
      KEGG.hyperG.matrix=summary(KEGG.hyperG.results)
      outMatrixName="kegg.hyperG.summary.csv"
      KEGG.hyperG.matrix$KEGGID=paste0('hsa:',KEGG.hyperG.matrix$KEGGID)
      KEGG.hyperG.matrix$p.adjust = p.adjust(KEGG.hyperG.matrix$Pvalue,method = 'BH')
      write.csv(KEGG.hyperG.matrix,outMatrixName )
    }

  }


}


