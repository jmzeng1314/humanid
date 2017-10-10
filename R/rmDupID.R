#' remove duplicat gene entrez ID or symbol for a expression matrix.
#'
#' sometimes, we will get a expression matrix, which one gene has many values in each sample, because we design many probeset for this  gene
#' most of time, we just need one value for each gene, so we should remove the duplicate value by choose the biggest.
#' the duplicted gene can be entrez gene ID or HUGO gene symbol
#'
#' @param dup_exprSet  a expression matrix or data.frame which the  first column is the ID needed to remove duplicate.
#' @return just  expression matrix which colname is the unique ID, and value is numeric.
#' @export
#' @keywords rmDupID
#' @examples
#' #' exprSet <- rmDupID(dup_exprSet)
#'
rmDupID <- function(dup_exprSet) {
    
    print(paste("Input is", nrow(dup_exprSet), "rows and after rmdup, just", length(unique(dup_exprSet[, 1])), "rows left", 
        sep = " "))
    
    exprSet = dup_exprSet[, -1]  ## first column is the ID needed to remove duplicate.probably is gene symbol or entrez ID.
    rowMeans = apply(exprSet, 1, function(x) mean(as.numeric(x), na.rm = T))
    dup_exprSet = dup_exprSet[order(rowMeans, decreasing = T), ]
    exprSet = dup_exprSet[!duplicated(dup_exprSet[, 1]), ]
    # exprSet=apply(exprSet,2,as.numeric)
    exprSet = exprSet[!is.na(exprSet[, 1]), ]
    rownames(exprSet) = exprSet[, 1]
    exprSet = exprSet[, -1]
    # str(exprSet)
    rn = rownames(exprSet)
    exprSet = apply(exprSet, 2, as.numeric)
    rownames(exprSet) = rn
    # exprSet[1:4,1:4]
    return(exprSet)
}
