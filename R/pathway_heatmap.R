#' draw heatmap for a single pathway or a batch of pahways
#'
#' you can point the pathway you are interested in or a list of related pathway ,such immune or sinal pathway
#'
#'
#' @param exprSet  a expression matrix, which columns are sample,rows are HUGO gene symbols.
#' @param pathwayID a character just like 00910(means 'Nitrogen metabolism'),default:00910
#' @param search  Generate heatmap for the kegg pathway which related the character you search.
#' @param all generate heatmap for all kegg pathway if set this parameter as TURE. default:F
#' @return generate QC figures for  expression matrix
#' @export
#' @keywords QC
#' @examples
#' pathway_heatmap(example_exprSet)


pathway_heatmap <- function(exprSet, pathwayID = "00910", search = "NA", all = F) {
    if (all) {
        drawPath = kegg2name
    } else if (search != "NA") {
        drawPath = kegg2name[grepl("signal", kegg2name$pathway_name), ]
    } else {
        if (is.na(match(pathwayID, kegg2name$pathway_id))) 
            stop("please input a correct pathway ID ,just like:00910")
        drawPath = kegg2name[pathwayID, ]
    }
    
    return("hello")
}
