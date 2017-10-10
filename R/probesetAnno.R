#' annotate entrez gene id for the given probeset
#'
#' somethins
#'
#' @param probeList a character vector for probeset ID list,default:c('1000_at','1001_at')
#' @param platform the platform for probeList, So far, we just accept hgu95av2/hgu133a/hgu133b/hgu133plus2
#' @return a data.frame which has two column, first is probeset ID , second is entrez gene id .
#' @export
#' @keywords microarray


probesetAnno <- function(probeList = c("1000_at", "1001_at"), platform = "hgu95av2") {
    if (length(probeList) > length(unique(probeList))) 
        warning("there is duplicate for the probeList !")
    if (platform == "hgu95av2") {
        results = hgu95av2_id[match(c(probeList), hgu95av2_id$probe_id), ]
    } else if (platform == "hgu133a") {
        results = hgu133a[match(c(probeList), hgu133a$probe_id), ]
    } else if (platform == "hgu133b") {
        results = hgu133b[match(c(probeList), hgu133b$probe_id), ]
    } else if (platform == "hgu133plus2") {
        results = hgu133plus2[match(c(probeList), hgu133plus2$probe_id), ]
    } else {
        stop("we do not accept this kind of microarry platform so far !")
    }
    return(results)
}






