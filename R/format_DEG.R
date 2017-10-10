#' write some files from the DEG results from topTable(limma)
#'
#'
#' @param DEG  DEG=topTable(fit,coef=2,adjust='BH')
#' @param prefix The prefix for all of the output files.
#' @param logFC_cutoff   The cutoff for logFC    to choose DEG, default:mean(abs(DEG$logFC)) + 2*sd(abs(DEG$logFC))
#' @param pvalue_cutoff  The cutoff for pvalue   to choose DEG, default:0.05
#' @param padjust_cutoff The cutoff for p.adjust to choose DEG, default:none
#' @return write some files
#' @export
#' @keywords format_DEG
#' @examples
#' #'  format_DEG(DEG)
#'
#'
format_DEG <- function(DEG, prefix = "test", logFC_cutoff = 0, pvalue_cutoff = 0, padjust_cutoff = 0, GOstats = F) {
    if (logFC_cutoff == 0) {
        logFC_cutoff <- mean(abs(DEG$logFC)) + 2 * sd(abs(DEG$logFC))
    }
    
    if (pvalue_cutoff == 0) {
        pvalue_cutoff <- 0.05
    }
    print(paste0("The cutoff for the logFC is : ", logFC_cutoff))
    print(table(abs(DEG$logFC) > logFC_cutoff & DEG$P.Value < pvalue_cutoff))
    
    if (padjust_cutoff == 0) {
        
        DEG$sigORnot <- ifelse(abs(DEG$logFC) > logFC_cutoff & DEG$P.Value < pvalue_cutoff, ifelse(DEG$logFC > logFC_cutoff, 
            "UP", "DOWN"), "NOT")
        
    } else {
        DEG$sigORnot <- ifelse(abs(DEG$logFC) > logFC_cutoff & DEG$P.Value < padjust_cutoff, ifelse(DEG$logFC > 
            logFC_cutoff, "UP", "DOWN"), "NOT")
        
    }
    print(table(DEG$sigORnot))
    
    Volcanic_DEG(DEG, imageType = "png", prefix = prefix)
    
    # png( paste0(prefix,'_volcanic.png') ) plot(DEG$logFC,DEG$P.Value,main =logFC_cutoff ) abline(h =
    # 0.05,col='blue') abline(v=logFC_cutoff,col='red') abline(v=-logFC_cutoff,col='red') dev.off()
    
    file_allGeneList = paste0(prefix, "_allGeneList.txt")
    file_diffGeneList = paste0(prefix, "_diffGeneList.txt")
    file_upGeneList = paste0(prefix, "_upGeneList.txt")
    file_downGeneList = paste0(prefix, "_downGeneList.txt")
    
    write.table(unique(DEG$symbol), file_allGeneList, quote = F, row.names = F, col.names = F)
    write.table(unique(DEG[DEG$sigORnot != "NOT", "symbol"]), file_diffGeneList, quote = F, row.names = F, col.names = F)
    write.table(unique(DEG[DEG$sigORnot == "UP", "symbol"]), file_upGeneList, quote = F, row.names = F, col.names = F)
    write.table(unique(DEG[DEG$sigORnot == "DOWN", "symbol"]), file_downGeneList, quote = F, row.names = F, col.names = F)
    
    batch_enrichment(file_diffGeneList, file_allGeneList, prefix = paste0(prefix, "_diff"), GOstats = GOstats)
    batch_enrichment(file_upGeneList, file_allGeneList, prefix = paste0(prefix, "_UP"), GOstats = GOstats)
    batch_enrichment(file_downGeneList, file_allGeneList, prefix = paste0(prefix, "_DOWN"), GOstats = GOstats)
    
    
    logFC <- DEG$logFC
    names(logFC) <- DEG$symbol
    
    keggEnrichTable = read.csv(paste0(prefix, "_diff_update_kegg.enrichment.csv"), stringsAsFactors = F)
    keggEnrichTable = keggEnrichTable[keggEnrichTable$Pvalue < 0.05, ]
    diff_gene_list = read.table(file_diffGeneList, stringsAsFactors = F)[, 1]
    add_kegg_up_down_link(keggEnrichTable, logFC, diff_gene_list, prefix = paste0(prefix, "_diff"))
    
    keggEnrichTable = read.csv(paste0(prefix, "_UP_update_kegg.enrichment.csv"), stringsAsFactors = F)
    keggEnrichTable = keggEnrichTable[keggEnrichTable$Pvalue < 0.05, ]
    diff_gene_list = read.table(file_diffGeneList, stringsAsFactors = F)[, 1]
    add_kegg_up_down_link(keggEnrichTable, logFC, diff_gene_list, prefix = paste0(prefix, "_UP"))
    
    keggEnrichTable = read.csv(paste0(prefix, "_DOWN_update_kegg.enrichment.csv"), stringsAsFactors = F)
    keggEnrichTable = keggEnrichTable[keggEnrichTable$Pvalue < 0.05, ]
    diff_gene_list = read.table(file_diffGeneList, stringsAsFactors = F)[, 1]
    add_kegg_up_down_link(keggEnrichTable, logFC, diff_gene_list, prefix = paste0(prefix, "_DOWN"))
    
    
    
    tmpAnno <- geneAnno(unique(DEG$symbol))
    dim(tmpAnno)
    DEG <- merge(DEG, tmpAnno, by = "symbol")
    dim(DEG)
    write.csv(DEG, paste0(prefix, "_DEG.csv"), row.names = F)
    
    
    rnk_file = paste0(prefix, ".rnk")
    sink(rnk_file)
    cat("#")
    sink()
    write.table(DEG[, c("symbol", "logFC")], rnk_file, append = T, quote = F, row.names = F)
    
    
}
