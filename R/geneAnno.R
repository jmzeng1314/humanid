#' annotate gene name or map for entrez id or symbol
#'
#' somethins
#'
#' @param geneLists a character vector for gene entrez ID list,default: c(1,2,9)
#' @param name     choose whether annotate gene name     for the gene id or not , default : T
#' @param map      choose whether annotate gene map      for the gene id or not , default : F
#' @param ensembl  choose whether annotate gene ensembl  for the gene id or not , default : F
#' @param accnum   choose whether annotate gene accnum   for the gene id or not , default : F
#' @return a data.frame which has more than 2 column
#' @export
#' @keywords geneAnno


geneAnno <- function(geneLists=c(1,2,9),name=T,map=F,ensembl=F,accnum=F){
  if(length(geneLists)>length(unique(geneLists)))
    warning('there is duplicate for the geneLists !')

  if( all(! geneLists %in% all_EG) ){
    inputType='symbol'
    geneLists=data.frame(symbol=geneLists)
    results=merge(geneLists,EG2Symbol,by='symbol',all.x=T)
  }else{
    inputType='entrezID'
    geneLists=data.frame(gene_id=geneLists)
    results=merge(geneLists,EG2Symbol,by='gene_id',all.x=T)
  }

  if ( name ){

    results=merge(results,EG2name,by='gene_id',all.x=T)
  }
  if(map){

    results=merge(results,EG2MAP,by='gene_id',all.x=T)
  }
  if(ensembl){

    results=merge(results,EG2ENSEMBL,by='gene_id',all.x=T)
  }
  if(accnum){

    results=merge(results,EG2MAP,by='gene_id',all.x=T)
  }
  results$symbol <- as.character(results$symbol  )
  return(results)

}





