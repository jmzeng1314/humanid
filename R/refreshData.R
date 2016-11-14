#' extract probeset ID and entrez gene id from R packages
#'
#' So far, we just extract data from R packages:hgu95av2/hgu133a/hgu133b/hgu133plus2
#'
#' @param refresh Don't use this function if not neccessary
#' @keywords refresh

refreshMicroarrayData <- function(refresh=F){
  if(refresh){
    library(hgu95av2.db)
    library(hgu133a.db)
    library(hgu133b.db)
    library(hgu133plus2.db)

    hgu95av2_id <- toTable(hgu95av2ENTREZID)
    hgu133a_id  <- toTable(hgu133aENTREZID)
    hgu133b_id  <- toTable(hgu133bENTREZID)
    hgu133plus2_id <- toTable(hgu133plus2ENTREZID)

    #save(hgu95av2_id,hgu133a_id,hgu133b_id,hgu133plus2_id,file = 'probeset.Rdata')
    devtools::use_data(hgu95av2_id,hgu133a_id,hgu133b_id,hgu133plus2_id, overwrite =T)
  }


  #hgu95av2_id[ match(c(probeList),hgu95av2_id$probe_id),]
}

#' extract gene IDs from org.Hs.eg.db
#'
#' So far, we just extract gene entrez ID,ensembl ID,HUGO symbol,and locus.
#'
#' @param refresh Don't use this function if not neccessary
#' @keywords refresh

refreshOrgHsData <- function(refresh=F){
  if(refresh){
    library(org.Hs.eg.db)
    all_EG=mappedkeys(org.Hs.egSYMBOL)
    allSymbols=mappedkeys(org.Hs.egSYMBOL2EG)
    EG2Symbol=toTable(org.Hs.egSYMBOL)
    EG2name=toTable(org.Hs.egGENENAME)
    EG2MAP=toTable(org.Hs.egMAP)
    EG2ENSEMBL=toTable(org.Hs.egENSEMBL)
    EG2accnum=toTable(org.Hs.egREFSEQ)
    devtools::use_data(all_EG,allSymbols,EG2Symbol,EG2name,EG2MAP,EG2ENSEMBL,EG2accnum, overwrite =T)
  }
}

#' update kegg data
#'
#' It's not a automately function, please don't run it. it will read the files ind data folder, which is 'kegg2geneID.txt' and 'kegg_hierarchical.txt'
#' you must make sure the files exists in the data/KEGG_update folder.
#' It will process the files and store the GeneID2kegg_list,kegg2GeneID_list,keggID2geneID_df,kegg2name to the Rdata.
#' GeneID2kegg_list and kegg2GeneID_list are list which one key has multiple values.
#' keggID2geneID_df is a data.frame , 2 columns : 'gene_id','path_id'
#' kegg2name is a data.frame, which has 4 columns, which is :parent1 parent2 pathway_id pathway_name
#'
#' @param refresh Don't use this function if not neccessary
#' @keywords refresh

update_kegg <- function(refresh=F){
  if(refresh){
    dir_bin='data'
    path2gene_file=file.path(dir_bin,'KEGG_update','kegg2geneID.txt')
    path2name_file=file.path(dir_bin,'KEGG_update','kegg_hierarchical.txt')
    if (file.exists(path2gene_file)){
      keggID2geneID=read.table(path2gene_file,sep="\t",colClasses=c('character'))
      keggID2geneID_df=keggID2geneID[,c(2,1)]
      names(keggID2geneID_df)=c('gene_id','path_id')
      tmp=read.table(path2gene_file,sep="\t",colClasses=c('character'))
      #tmp=toTable(org.Hs.egPATH)
      # first column is kegg ID, second column is entrez ID
      GeneID2kegg_list<<- tapply(tmp[,1],as.factor(tmp[,2]),function(x) x)
      kegg2GeneID_list<<- tapply(tmp[,2],as.factor(tmp[,1]),function(x) x)
    }else{stop("we can not find the file:path2gene_file")}
    if (file.exists(path2name_file)){
      kegg2name<<- read.delim(path2name_file,header=F,sep="\t",colClasses=c('character'),stringsAsFactors =F)
      colnames(kegg2name)=c('parent1','parent2','pathway_id','pathway_name')
      ###kegg2name$pathway_id=as.numeric(kegg2name$pathway_id)
      rownames(kegg2name)=kegg2name$pathway_id
    }else{stop("we can not find the file:path2name_file")}

    devtools::use_data(keggID2geneID_df,GeneID2kegg_list,kegg2GeneID_list, kegg2name,overwrite =T)
  }
}

#' create the duplicated expression matrix
#'
#' the expression matrix come from the data 'sCLLex' in R package 'CLL'
#' duplicated because one gene symbol corresponding to multiple probeset ID in hgu95av2.db
#' dup_exprSet is a data.frame,  the first column is gene symbol, which needed to remove duplicated.
#'
#' @param refresh Don't use this function if not neccessary
#' @keywords refresh
#'
create_dup_exprSet <- function(refresh=F){
  library(CLL)
  data(sCLLex)
  group_list=sCLLex$Disease

  library(annotate)
  exprSet=exprs(sCLLex)
  platformDB='hgu95av2.db'
  library(platformDB, character.only=TRUE)
  probeset <- featureNames(sCLLex)
  SYMBOL <-  lookUp(probeset, platformDB, "SYMBOL")
  dup_exprSet=cbind(SYMBOL,exprSet)

  devtools::use_data(dup_exprSet,overwrite =T)
}

#' create expression matrix,rowname are unique genes, colname are samples
#'
#' example_exprSet comes from rmDupID(dup_exprSet),also the data 'sCLLex' in R package 'CLL'
#'
#' @param refresh Don't use this function if not neccessary
#' @keywords refresh
#'
create_example_exprSet <- function(refresh=F){
  example_exprSet=rmDupID(dup_exprSet)
  devtools::use_data(example_exprSet,overwrite =T)
}


#' create a list of gene sets for  enzyme/protein complex/ development/immune
#'
#' PTM: Post-translational modification, such methyl-,  acethyl-,  phospha-
#' SWI/SNF Chromatin Remodeling Complexes and PRC family members
#' development associated genes
#' Immune-Related Genes
#' transcription factor genesets
#' other genesets
#'
#' @param refresh Don't use this function if not neccessary
#' @keywords refresh
#'
create_genesets_list <- function(refresh=F){

  suppressMessages(library("org.Hs.eg.db"))
  allSymbols=mappedkeys(org.Hs.egSYMBOL2EG)
  PRMT_list =allSymbols[grepl("^PRMT", allSymbols)];PRMT_list=c(PRMT_list,'CARM1')
  JMJD_list =allSymbols[grepl("^JMJD", allSymbols)];
  KMT_list  =allSymbols[grepl("^KMT",  allSymbols)];KMT_list=c(KMT_list,'EZH2')
  KDM_list  =allSymbols[grepl("^KDM",  allSymbols)];
  HAT_list  =allSymbols[grepl("^KAT\\d",  allSymbols)];HAT_list=c(HAT_list,'HAT1','EP300')
  HDAC_list =allSymbols[grepl("^HDAC",  allSymbols)];
  MAPK_list =allSymbols[grepl("^MAPK",  allSymbols)];
  CDK_list  =allSymbols[grepl("^CDK\\d+$", perl =T, allSymbols)];
  PRAK_list  =allSymbols[grepl("^PRKA", perl =T, allSymbols)];
  DUSP_list =allSymbols[grepl("^DUSP\\d+$",  allSymbols)];

  enzyme_genesets=list(PRMT_list,JMJD_list,KMT_list,KDM_list,HAT_list,HDAC_list,MAPK_list,CDK_list,PRAK_list,DUSP_list)
  names(enzyme_genesets)=strsplit('PRMT_list,JMJD_list,KMT_list,KDM_list,HAT_list,HDAC_list,MAPK_list,CDK_list,PRAK_list,DUSP_list',',')[[1]]

  HIST_list =allSymbols[grepl("^HIST",  allSymbols)];
  ABC_list =allSymbols[grepl("^ABC",  allSymbols)];
  TUB_list =allSymbols[grepl("^TUB",  allSymbols)];
  ACT_list =allSymbols[grepl("^ACT",  allSymbols)];
  MYO_list =allSymbols[grepl("^MYO",  allSymbols)];
  other_genesets=list(HIST_list,ABC_list,TUB_list,ACT_list,MYO_list)
  names(other_genesets)=strsplit('HIST_list,ABC_list,TUB_list,ACT_list,MYO_list',',')[[1]]


  PRC_list = c('BMI1','CBX2','CBX4','CBX7','CBX8','PCGF1','PCGF2','PCGF3','PCGF5','PCGF6','PHC1','PHC2','PHC3','RING1','RNF2','SCML1','SCML2','AEBP2','EED','EZH2','HOXTAIR','JARID2','RBBP4','RBBP6','SUZ12','GAPDH')
  SWI_SNF_list= allSymbols[grepl("^SMAR",  allSymbols)];
  RNA_Polymerase_II_list = allSymbols[grepl("^POLR",  allSymbols)];
  TAF_list =allSymbols[grepl("^TAF",  allSymbols)];
  TBP_list =allSymbols[grepl("^TBP",  allSymbols)];
  TATA_Box_TF_list=c(TAF_list,TBP_list)
  General_TF_list =allSymbols[grepl("^GTF",  allSymbols)];
  nulear_Factor_list =allSymbols[grepl("^NF",  allSymbols)];
  protein_complex_genesets=list(PRC_list,SWI_SNF_list,RNA_Polymerase_II_list,TATA_Box_TF_list,General_TF_list,nulear_Factor_list)
  names(protein_complex_genesets)=strsplit('PRC_list,SWI_SNF_list,RNA_Polymerase_II_list,TATA_Box_TF_list,General_TF_list,nulear_Factor_list',',')[[1]]


  devtools::use_data(protein_complex_genesets,enzyme_genesets,other_genesets,overwrite =T)

}







