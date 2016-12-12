#' Do basic QC for a expression matrix.
#'
#' Basic QC include:boxplot and hist figures for expression value distribution
#'                  correlations among the samples by clustering
#'
#'
#' @param exprSet  a expression matrix, which columns are sample,rows are HUGO gene symbols.
#' @param group_list a vector,as long as the col number for the expression matrix,which describe the group for the samples in exprSet
#' @param prefix The prefix for all of the output files.
#' @return generate QC figures for  expression matrix
#' @export
#' @keywords QC
#' @examples
#' #' QCexpressionMatrix(example_exprSet,group_list,prefix='test')



QCexpressionMatrix <- function(exprSet,group_list='NO',prefix='test'){
  library(reshape2)
  library(ggplot2)

  #xpr_lng=reshape2::melt(exprSet)
  exprSet_L=melt(exprSet)
  colnames(exprSet_L)=c('probe','sample','value')
  p=ggplot(exprSet_L,aes(x=sample,y=value))+geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
  png(paste0(prefix,"_values_distribution_boxplot.png"))
  print(p)
  dev.off()

  png(paste0(prefix,"_values_distribution_hist.png"))
  p=ggplot(exprSet_L,aes(colour =sample,x=value))+geom_density()
  print(p)
  dev.off()

  png(paste0(prefix,"_sample_correlation.png"))
  out.dist=dist(t(exprSet),method='euclidean')
  out.hclust=hclust(out.dist,method='complete')
  plot(out.hclust)
  dev.off()

  if ( length(group_list) >1){

    xpr=exprSet
    this.color=rainbow(length(unique(group_list)))

    ## bug:cannot rescale a constant/zero column to unit variance
    ## https://searchcode.com/codesearch/view/15421003/
    xpr_t = t(xpr)
    variances <- apply(xpr_t, 2, var) # variance for each column
    ## some probese/gene show same expression values in all of the samples
    zerovar <- which(variances == 0)
    if (length(zerovar)==0){
      xpr_t.2 <- xpr_t
    }else{

      xpr_t.2 <- xpr_t[,-zerovar] # delete columns where variance is zero
    }
    ## column mean gene and row means sample
    pc <- prcomp(xpr_t.2,scale=TRUE)
    ##pc=prcomp(t(xpr), scale=T)


    #summary(pc)
    #screeplot(pc, type='line')
    pcx=data.frame(pc$x)
    pcr=cbind(samples=rownames(pcx),group_list, pcx)
    png(paste0(prefix,"_sample_PCA.png"),width = 800,height = 800)
    p=ggplot(pcr, aes(PC1, PC2))+geom_point(size=5, aes(color=group_list))+scale_colour_manual(values=this.color)+
      geom_text(aes(label=samples),hjust=-0.1, vjust=-0.3)
    print(p)
    dev.off()


    xpr=exprSet
    this.color=rainbow(length(unique(group_list)))
    #as.numeric(factor(group_list))
    meta=data.frame(Group=group_list,
                    Sample=colnames(exprSet),
                    Color=this.color[as.numeric(factor(group_list))],
                    stringsAsFactors = F
    )
    tmp=cc_plot(xpr, meta)
    png(paste0(prefix,"_sample_cc.png"),width = 800,height = 800)
    print(tmp$plt)
    dev.off()
    write.csv(tmp$dfp,paste0(prefix,"_sample_cc.csv") ,quote = F)


  }
}



