#' annotate KEGG pathway for entrez id or symbol
#'
#' somethins
#'
#' @param geneLists a character vector for gene entrez ID list,default: c(1,2,9)
#' @param fileType     choose whether annotate gene name     for the gene id or not , default : T
#' @param fileName      choose whether annotate gene map      for the gene id or not , default : F
#' @param multiple  choose whether annotate gene ensembl  for the gene id or not , default : F
#' @return write a html or csv files for the annotation results.
#' @export
#' @keywords keggAnno
#' @examples
#' #' keggAnno('TP53',fileName='TP53',multiple = F),keggAnno('TP53')



keggAnno <- function(geneLists=c(1,2,9),fileType='html',fileName='gene',multiple=T){
  results <- geneAnno(geneLists)
  tmp=merge(results,keggID2geneID_df,by = 'gene_id', all.x=T)
  keggID2name=kegg2name[,3:4]
  colnames(keggID2name)=c('path_id','path_name')
  tmp2=merge(tmp,keggID2name,by='path_id',all.x=T)
  kegg_info=tmp2
  ## 5 columns: path_id,gene_id,symbol,gene_name,path_name
  ## we need to add links for the path according to the path_id and path_name !!

  if(fileType !='html'){
    write.csv(kegg_info, paste0(fileName,'_kegg_anno.csv'),row.names = F)
  }
  else{
    if(multiple){
      kegg_info$pathway_name=apply(kegg_info,1,function(x) {
        ## create a new column for kegg_info !!!
        ##path_id gene_id symbol path_name
        tmp='';
        if (is.na(x['path_id'])){
          tmp='NA'
        }else{

          pathwayID=sprintf("%05.0f",as.numeric(x['path_id']))
          pathwayName=x['path_name']
          href=paste0("http://www.genome.jp/kegg-bin/show_pathway?hsa",
                      pathwayID,"+hsa:",x['gene_id']
          )
          tmp=paste0('<b><a target="_black" href=', shQuote(href) ,'>',pathwayName,'</a></b>')

        }
        return(tmp)
      }) ## we add the links in the new column: pathway_name
      new_kegg_info=kegg_info[order(kegg_info[,'gene_id']),]
    }else{
      new_kegg_info <- data.frame(matrix(unlist(lapply(split(kegg_info,kegg_info$gene_id), function(gene_id_kegg_info){
        gene_id_kegg_info$pathway_name <- apply(gene_id_kegg_info,1,function(x) {
          tmp='';
          if (is.na(x['path_id'])){
            tmp='NA'
          }else{

            pathwayID=sprintf("%05.0f",as.numeric(x['path_id']))
            pathwayName=x['path_name']
            href=paste0("http://www.genome.jp/kegg-bin/show_pathway?hsa",
                        pathwayID,"+hsa:",x['gene_id']
            )
            tmp=paste0('<b><a target="_black" href=', shQuote(href) ,'>',pathwayName,'</a></b>')

          }
          return(tmp)
        }) ## end for apply
        links=paste(gene_id_kegg_info$pathway_name,collapse=',')
        return(c(gene_id_kegg_info[1,c("gene_id","symbol","gene_name")],links))
      })## end for lapply +split
      ) ## end for unlist
      , ncol=4, byrow=T) ## end for matrix
      ,stringsAsFactors=FALSE) ## endor data.frame
      colnames(new_kegg_info)=c("gene_id","symbol","gene_name",'kegg_links')
    }
    y <- DT::datatable(new_kegg_info,escape = F,rownames=F)
    DT::saveWidget(y, paste0(fileName,'_kegg_anno.html'))
  }

}
