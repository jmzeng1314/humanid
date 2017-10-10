#' annotate gene name or map for entrez id or symbol
#'
#' So far, we accept entrez gene id or symbol as input to be annotated. such as 7157 or TP53. you just need to assign the geneLists,a character vector
#' Then we'll annotate the name by default, and also you can choose map,ensembl,accnum to be annotated or not.
#' lastly, you can output the annotation to a html file,just like:geneAnno(kegg2GeneID[['01212']],file = T,prefix = as.character(kegg2name['01212','pathway_name'])
#' for each kegg pathway, you will like this function I guess. So have a fun .
#'
#' @param geneLists a character vector for gene entrez ID list,default: c(1,2,9)
#' @param name     choose whether annotate gene name     for the gene id or not , default : T
#' @param map      choose whether annotate gene map      for the gene id or not , default : F
#' @param ensembl  choose whether annotate gene ensembl  for the gene id or not , default : F
#' @param accnum   choose whether annotate gene accnum   for the gene id or not , default : F
#' @param file     choose whether print the annotate results into a file or not , default : F
#' @param prefix   define the prefix for the output file, default : test
#' @return a data.frame which has more than 2 column
#' @export
#' @keywords geneAnno
#' @examples
#' geneAnno();geneAnno('TP53');geneAnno(kegg2GeneID[['01212']],file = T,prefix = as.character(kegg2name['01212','pathway_name'])
#' lapply(names(kegg2GeneID),function(x) geneAnno(kegg2GeneID[[x]],file = T,prefix = as.character(kegg2name[x,'pathway_name'])))
#'

geneAnno <- function(geneLists = c(1, 2, 9), name = T, map = F, ensembl = F, accnum = F, file = F, prefix = "test") {
    if (length(geneLists) > length(unique(geneLists))) 
        warning("there is duplicate for the geneLists !")
    
    geneLists <- unique(geneLists)
    if (all(!geneLists %in% all_EG)) {
        inputType = "symbol"
        geneLists = data.frame(symbol = geneLists)
        results = merge(geneLists, EG2Symbol, by = "symbol", all.x = T)
    } else {
        inputType = "entrezID"
        geneLists = data.frame(gene_id = geneLists)
        results = merge(geneLists, EG2Symbol, by = "gene_id", all.x = T)
    }
    
    if (name) {
        
        results = merge(results, EG2name, by = "gene_id", all.x = T)
    }
    if (map) {
        
        results = merge(results, EG2MAP, by = "gene_id", all.x = T)
    }
    if (ensembl) {
        
        results = merge(results, EG2ENSEMBL, by = "gene_id", all.x = T)
    }
    if (accnum) {
        
        results = merge(results, EG2MAP, by = "gene_id", all.x = T)
    }
    results$symbol <- as.character(results$symbol)
    if (file) {
        entrez_prefix <- "http://www.ncbi.nlm.nih.gov/gene/"
        href = paste0(entrez_prefix, results$gene_id)
        results$gene_id = paste0("<b><a target=\"_black\" href=", shQuote(href), ">", results$gene_id, "</a></b>")
        
        symbol_prefix <- "http://www.ncbi.nlm.nih.gov/gene?term="
        href = paste0(symbol_prefix, results$symbol)
        results$symbol = paste0("<b><a target=\"_black\" href=", shQuote(href), ">", results$symbol, "</a></b>")
        # prefix <- 'fefe / \ ' fee'
        prefix <- sub("/", "_", prefix)
        prefix <- sub("\"", "_", prefix)
        
        y <- DT::datatable(results, escape = F, rownames = F)
        DT::saveWidget(y, paste0(prefix, "_geneAno_links.html"))
    } else {
        return(results)
    }
    
    
}





