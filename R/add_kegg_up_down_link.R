#' add links for the KEGG enrichment table.
#'
#'  KEGG enrichment table comes from the stantard hyperG test. the PathwayID and pathway_name columns  are needed.
#'  logFC <- DEG$logFC;names(logFC) <- DEG$symbol
#'  diff_gene_list also comes from DEG results.
#'
#' @param keggEnrichTable a data.frame which must have PathwayID and pathway_name columns
#' @param logFC  A numeric vector, the names for the vector are gene lists.
#' @param diff_gene_list   a vector which contain the significantly DEGs. symbol or entrez ID.
#' @param prefix The prefix for all of the output files.
#' @return write a html file with links for the keggEnrichTable
#' @export
#' @keywords keggAnno
#' @examples
#' #' add_kegg_up_down_link(keggEnrichTable , DEG , diff_gene_list,prefix)



add_kegg_up_down_link <- function(keggEnrichTable, logFC, diff_gene_list, prefix = "test") {
    keggEnrichTable$links = apply(keggEnrichTable, 1, function(x) {
        this_keggID <- sub("hsa:", "", x["PathwayID"])
        
        this_kegg_has_geneID <- kegg2GeneID_list[[this_keggID]]
        this_kegg_has_geneSymbol <- unique(geneAnno(this_kegg_has_geneID)$symbol)
        this_kegg_has_geneSymbol_diff <- intersect(this_kegg_has_geneSymbol, diff_gene_list)
        ## UP GENE: DHX8+red%2Cblue%0D DOWN GENE:DHX8+blue%2Cred%0D
        color = ifelse(logFC[this_kegg_has_geneSymbol_diff] > 0, "red%2Cblue%0D", "blue%2Cred%0D")
        kegg_suffix <- paste0(paste(this_kegg_has_geneSymbol_diff, color, sep = "+"), collapse = "")
        
        href = paste0("http://www.genome.jp/kegg-bin/show_pathway?map=hsa", this_keggID, "&multi_query=", kegg_suffix)
        link = paste0("<b><a target=\"_black\" href=", shQuote(href), ">", x["pathway_name"], "</a></b>")
        return(link)
        
    })
    
    y <- DT::datatable(keggEnrichTable, escape = F, rownames = F)
    DT::saveWidget(y, paste0(prefix, "_kegg_links.html"))
    
}
