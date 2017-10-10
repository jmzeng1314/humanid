#' create gct file and cls files
#'
#' create gct file and cls file for GSEA according the expression matrix and group information
#'
#' @param prefix The prefix for gct file and cls files.
#' @param exprSet a expression matrix, which columns are sample,rows are HUGO gene symbols, or probeset ID .
#' @param group_list a vector,as long as the col number for the expression matrix,which describe the group for the samples in exprSet
#' @param destdir where to store the files just download.
#' @return write 2 files which are the input for GSEA (gct and cls format)
#' @export
#' @keywords downGSE
#' @examples
#' #' createGSEAinput('GSE1009',exprSet,group_list)


createGSEAinput <- function(prefix = "GSE1009", exprSet = example_exprSet, group_list, destdir = ".") {
    # sink('outfile.txt') cat('hello') cat('\n') cat('world') sink()
    gct_file = paste0(prefix, ".gct")
    sink(gct_file)
    cat("#1.2\n")
    cat(paste0(nrow(exprSet), "\t", length(group_list), "\n"))
    sink()
    gct_out <- cbind(symbol = rownames(exprSet), description = "na", exprSet)
    write.table(gct_out, gct_file, append = T, quote = F, row.names = F, sep = "\t")
    
    
    cls_file = paste0(prefix, ".cls")
    sink(cls_file)
    cat(paste0(length(group_list), " ", length(unique(group_list)), " 1\n"))
    cat(paste0("# ", paste(unique(group_list), collapse = " "), "\n"))
    cat(paste(group_list, collapse = " "))
    sink()
}



