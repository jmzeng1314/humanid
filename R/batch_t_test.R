#' do batch t.test for a matrix
#'
#'
#' @param dat An expression matrix which rownames are genes and columns are samples
#' @param group1 defined which columns belong to the first group as control
#' @param group2 defined which columns belong to the second group
#' @return a data.frame just like topTable
#' @export
#' @keywords batch_enrichment
#' @examples
#' #'  batch_t_test(exprSet ,4:7,1:3)
#'
#'
batch_t_test <- function(dat=matrix(rnorm(300),30,10),group1=1:5,group2=6:10){

  dat1=dat[,group1];dat2=dat[,group2]
  dat=cbind(dat1,dat2)
  library(pi0)
  pvals=matrix.t.test(dat,1,length(group1),length(group2))
  p.adj=p.adjust(pvals,method = 'BH')
  avg_1=rowMeans(dat1);avg_2=rowMeans(dat2);
  FC=avg_2/avg_1;
  results=cbind(avg_1,avg_2,FC,pvals,p.adj)
  rownames(results)=rownames(dat)
  results=as.data.frame(results)
  return(results[order(results$pvals),])

}
