#' do batch t.test for a matrix
#'
#'

#' @param dat    Matrix with microarray expression values.
#' @param group_list Factors for two groups that are tested for differential expression.
#' @return a data.frame just like topTable
#' @export
#' @keywords batch_t_test
#' @examples
#' #'  batch_t_test()
#'
#'
batch_t_test <- function(dat=matrix(rnorm(300),30,10),group_list=factor(c(rep(1,5),rep(2,5))) ){

  group1= which(group_list==levels(group_list)[1])
  group2= which( group_list==levels(group_list)[2])
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
