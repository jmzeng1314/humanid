#' Hypergeometric Tests for GO/KEGG test by jimmy
#'
#' Given GeneID2Path,Path2GeneID,diff_gene,universeGeneIds, this function will compute  Hypergeometric Tests for each path (GO/KEGG),
#' It can't hold on the structure of the GO graph, just a simple path .
#'
#' @param GeneID2Path a list which one entrez gene id to multiple pathway id,default:GeneID2kegg_list
#' @param Path2GeneID a list which one pathway id to multiple entrezgene id,default:kegg2GeneID_list
#' @param diff_gene   a vector which contain the significantly DEG list, abouth 500~1000 genes(entrez gene id ),default:sample(unique(hgu95av2_id$gene_id),500)
#' @param universeGeneIds a vector which contain the backgroud gene list,probably 20,000 genes , default:unique(hgu95av2_id$gene_id)
#' @return a data.frame, each row is a Hypergeometric Tests result for each pathway .
#' @export
#' @keywords hyperGtest
#' @examples
#' #' hyperGtest_jimmy(),hyperGtest_jimmy()


hyperGtest_jimmy <- function(GeneID2Path = GeneID2kegg_list, Path2GeneID = kegg2GeneID_list, diff_gene = sample(unique(hgu95av2_id$gene_id), 
    500), universeGeneIds = unique(hgu95av2_id$gene_id)) {
    diff_gene_has_path = intersect(diff_gene, names(GeneID2Path))
    n = length(diff_gene)  #306
    N = length(universeGeneIds)  #5870
    results = c()
    
    for (i in names(Path2GeneID)) {
        M = length(intersect(Path2GeneID[[i]], universeGeneIds))
        # print(M)
        if (M < 5) 
            next
        exp_count = n * M/N
        # print(paste(n,N,M,sep='\t'))
        k = 0
        for (j in diff_gene_has_path) {
            if (i %in% GeneID2Path[[j]]) 
                k = k + 1
        }
        OddsRatio = k/exp_count
        if (k == 0) 
            next
        p = phyper(k - 1, M, N - M, n, lower.tail = F)
        # print(paste(i,p,OddsRatio,exp_count,k,M,sep=' '))
        results = rbind(results, c(i, p, OddsRatio, exp_count, k, M))
    }
    colnames(results) = c("PathwayID", "Pvalue", "OddsRatio", "ExpCount", "Count", "Size")
    results = as.data.frame(results, stringsAsFactors = F)
    results$p.adjust = p.adjust(results$Pvalue, method = "BH")
    results = results[order(results$Pvalue), ]
    rownames(results) = 1:nrow(results)
    return(results)
}
