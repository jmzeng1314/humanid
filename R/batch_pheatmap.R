#' Do batch heatmap for interested gene sets.
#'
#' extract the corresponding expression matrix for each gene set , and draw heatmap for the gene set expression matrix.
#' The name of the gene set will be used for the name of the heatmap figure
#' the rowname must be HUGO gene symbols
#'
#' @param exprSet  a expression matrix, which columns are sample,rows are HUGO gene symbols.
#' @param needGeneList a list of interested gene sets, each gene set is a vector of genes, the name will be used by out figures.
#' @return generate heatmap figures for each gene sets.
#' @export
#' @keywords heatmap
#' @examples
#' #' batch_pheatmap(),batch_pheatmap()


batch_pheatmap <- function(exprSet,needGeneList){
  fileList=names(needGeneList)
  lapply(1:length(fileList), function(i){
    x=needGeneList[[i]]
    x=x[ x %in% rownames(exprSet) ]
    print(x)
    if(length(x)>4){
      data1 = exprSet[x,]
      geneName=geneAnno(rownames(data1))
      geneName= geneName[match(x,geneName$symbol ),]
      geneName=paste0(geneName[,2],'(',geneName[,3],')')

      library(pheatmap)
      # color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
      if(length(x) > 40 ){
        pheatmap(data1,
                 cluster_rows=F,cluster_cols=T, #display_numbers = TRUE
                 border_color="white" ,labels_row = geneName, cellwidth=20,fontsize=6,
                 file=paste0('heatmap_',fileList[i],".pdf"),width=10,height=10
        )
      }else{
        png(paste0('heatmap_',fileList[i],".png"),width = 600,height = 600);
        pheatmap(data1,
                 cluster_rows=F,cluster_cols=T, #display_numbers = TRUE
                 border_color="white" ,labels_row = geneName
        )
      }
      dev.off()
    }
  })
}

