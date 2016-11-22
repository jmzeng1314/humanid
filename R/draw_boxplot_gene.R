#' draw boxplot for a specific gene's values between two groups
#'
#'
#' @param gene  specific a gene which should be the rownames of exprSet
#' @param exprSet An expression matrix which rownames are genes and columns are samples
#' @param group_list defined which sample belong to which group
#' @return create a boxplot
#' @export
#' @keywords batch_enrichment
#' @examples
#' #'  draw_boxplot_gene(x,exprSet,group_list)
#'
#'

draw_boxplot_gene <- function(gene ='VCX3A',exprSet,group_list){
  if (length(unique(group_list)) != 2)
    stop("we just accept 2 groups for the group_list!!!")
  if( gene  %in% rownames(exprSet)  ){
    png(paste0(gene,'_boxplot.png'))
    dat1= data.frame(value = as.numeric(exprSet[gene,]),
                     type = group_list
    )
    boxplot( value ~  type, data = dat1, lwd = 2, ylab = 'value')
    title(main = gene,sub = t.test(value ~  type, data = dat1)$p.value)
    stripchart(value ~ type, vertical = TRUE, data = dat1,
               method = "jitter", add = TRUE, pch = 20, col = 'blue')
    dev.off()
  }
}
# CTlist=read.table('CTlist.txt',stringsAsFactors = F)[,1]
# lapply(CTlist, function(x) draw_boxplot_gene(x,exprSet,group_list)  )

