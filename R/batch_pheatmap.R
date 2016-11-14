#' Do batch heatmap for interested gene sets.
#'
#' extract the corresponding expression matrix for each gene set , and draw heatmap for the gene set expression matrix.
#' The name of the gene set will be used for the name of the heatmap figure
#' the rowname must be HUGO gene symbols
#' we have prepared 3 genesets:protein_complex_genesets,enzyme_genesets,other_genesets
#'
#' @param exprSet  a expression matrix, which columns are sample,rows are HUGO gene symbols.
#' @param group_list a vector,as long as the col number for the expression matrix,which describe the group for the samples in exprSet
#' @param name Whether to attach the gene name to the heatmap or not, default:F
#' @param genesets a list of interested gene sets, each gene set is a vector of genes, the name will be used by out figures.
#' @return generate heatmap figures for each gene sets.
#' @export
#' @keywords heatmap
#' @examples
#' #' batch_pheatmap(),batch_pheatmap()


batch_pheatmap <- function(exprSet,group_list,name=F,genesets,width=800,height=800){
  fileList=names(genesets)
  sample_list=colnames(exprSet)
  lapply(1:length(fileList), function(i){
    x=genesets[[i]]
    x=x[ x %in% rownames(exprSet) ]
    print(x)
    if(length(x)>4){
      data1 = exprSet[x,]

      ## Need to be simpler,but I can't do it.
      tmp_number=rep(0,length(unique( group_list )))
      names(tmp_number)=unique( group_list )
      tmp_sample_name=c()
      for(j in 1:length( sample_list )){
        tmp_number[ group_list[j] ]= tmp_number[ group_list[j] ] +1;
        tmp_sample_name=c(tmp_sample_name,paste(group_list[j], tmp_number[ group_list[j] ]  ,sep="_"))
      }
      #print(tmp_sample_name)

      data1=data1[,sample_list[order(tmp_sample_name)]]

      geneName = rownames(data1)
      if(name){
        geneName=geneAnno(rownames(data1))
        geneName= geneName[match(x,geneName$symbol ),]
        geneName=paste0(geneName[,2],'(',geneName[,3],')')
      }

      # Generate annotations for rows and columns
      annotation_col = data.frame(
        sampleGroup = group_list
      )
      rownames(annotation_col) =  sample_list



      library(pheatmap)
      # color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
      if(length(x) > 40 ){
        pheatmap(data1,
                 cluster_rows=T,cluster_cols=T, annotation_col = annotation_col,
                 cutree_col = 2, #display_numbers = TRUE
                 border_color="white" ,labels_row = geneName, cellwidth=20,fontsize=6,
                 file=paste0('heatmap_',fileList[i],".pdf"),width=10,height=10
        )
      }else{
        png(paste0('heatmap_',fileList[i],".png"),width = width,height = height);
        pheatmap(data1,
                 cluster_rows=T,cluster_cols=T, annotation_col = annotation_col,
                 cutree_col = 2, #display_numbers = TRUE
                 border_color="white" ,labels_row = geneName
        )
      }
      dev.off()
    }
  })
}




