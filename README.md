# Jimmy's first R package:humanid

------

It's a unoffical package, just for **fun**. I will put many simple R functions which used frequently in my daily data analysis, which should be useful for you too. Don't be worry to learn them, I'm used to write minimal function, just like perl, my favorate computer language ! 

### functions including:

> * batch_enrichment
> * batch_pheatmap
> * batch_t_test
> * createGSEAinput
> * downGSE
> * draw_boxplot_gene
> * format_DEG
> * geneAnno
> * keggAnno
> * pathway_heatmap(todo)
> * probesetAnno
> * QCexpressionMatrix 


please feel free to contact with me if there's bug in my package, so far it's very simple，not Defensive programming,So don't fool around with my function. 

### [my blog](http://www.bio-info-trainee.com/)

>  I perfer email communication: jmzeng1314@163.com 

------

## A basic example to use this package:

```R
lapply(c('humanid','GEOquery','hugene10sttranscriptcluster.db'),function(x) library(x,character.only = T))
studyID <- 'GSE42872'
eSet<-getGEO(studyID, destdir='./',getGPL= F)
## you'd check the pdata(eSet) manually 
exprSet<-get_symbol_exprSet(eSet[[1]],'hugene10sttranscriptcluster.db')
QCexpressionMatrix(exprSet = exprSet,group_list =c(1,1,1,0,0,0) ,prefix =studyID )
createGSEAinput(prefix =studyID ,exprSet,group_list =c(1,1,1,0,0,0))
DEG <-do_DEG_2groups(prefix =studyID ,
                     exprSet= exprSet,
                     group_list=c(1,1,1,0,0,0),
                     method='limma')
DEG$symbol<-rownames(DEG);
format_DEG(DEG,prefix =studyID ,GOstats = T) 
## it will take about 10 minutes, you can set the GOstats = F to speed up this process.
```



```R
## this is a test, don't run it.
source("https://bioconductor.org/biocLite.R")
biocLite("zebrafishRNASeq")
library(zebrafishRNASeq)
data(zfGenes)
head(zfGenes)
spikes <- zfGenes[grep("^ERCC", rownames(zfGenes)),]
head(spikes)

## below is the example code.

library(CLL)
data(sCLLex)
group_list=sCLLex$Disease
QCexpressionMatrix(example_exprSet,group_list)
batch_pheatmap(example_exprSet,group_list,name=F,genesets=enzyme_genesets)

rm(list=ls())
library(humanid)
library(GEOquery)
library(affy)
library(limma)
old_wd='G:/array/GSE34824'
setwd(old_wd)
studyID='GSE34824';
studyID_probe=paste0(studyID,'_probe')
studyID_gene=paste0(studyID,'_gene')

R_history_data <- paste0(studyID,'.Rdata')
if ( file.exists(R_history_data)){
  load( R_history_data )
}else{
  eSet <-  downGSE(studyID)
  save.image( R_history_data )
}

exprSet=exprs(eSet[[1]])
pdata=pData(eSet[[1]])

treatment=factor(unlist(lapply( pdata$characteristics_ch1.3 ,function(x) strsplit(as.character(x),": ")[[1]][2]))) 
treatment=relevel(treatment,'wt')

## we just need to compare the K27M to the WT samples.

choose_eSet = eSet[[1]][,treatment %in% c('K27M','wt') ]
exprSet=exprs( choose_eSet )
pdata=pData( choose_eSet )

treatment=factor(unlist(lapply( pdata$characteristics_ch1.3 ,function(x) strsplit(as.character(x),": ")[[1]][2]))) 
treatment=relevel(treatment,'wt')

## create gct file and cls file for GSEA 
createGSEAinput(studyID , exprSet,treatment) 
group_list = treatment


## create a new folder to store the DEG results based on probeset ID 
setwd(old_wd)
if ( ! file.exists('basedonProbe') )
  dir.create('basedonProbe')

setwd('basedonProbe') 
studyID <- studyID_probe
if(T){  
  
  exprSet <- exprs(choose_eSet);dim(exprSet)
  exprSet <- na.omit(exprSet);dim(exprSet)
  #exprSet=as.data.frame(exprSet);dim(exprSet)
  exprSet <- na.omit(exprSet);dim(exprSet)
  if( mean(rowMeans( exprSet ,na.rm = T),na.rm = T) >20)
    exprSet=log2(exprSet) ## based on 2;
  dim(exprSet)
  
  
  QCexpressionMatrix(exprSet, group_list, project = studyID)
  
  library(limma)
  design=model.matrix(~ group_list )
  fit=lmFit(exprSet,design)
  fit=eBayes(fit)
  options(digits = 4)
  DEG <- topTable(fit,coef=2,adjust='BH',n=Inf) ;dim(DEG)
  
  DEG$probe_id = rownames(DEG)
  library(hgu133plus2.db)
  probe2symbol_df <- AnnotationDbi::toTable(hgu133plus2SYMBOL)
  DEG <- merge(probe2symbol_df,DEG,by='probe_id');dim(DEG)
  format_DEG(DEG,studyID)
  
}

setwd(old_wd)
if ( ! file.exists('basedonGene') )
  dir.create('basedonGene')

setwd('basedonGene')
studyID <- studyID_gene

exprSet <- get_symbol_exprSet(choose_eSet)
exprSet <- na.omit(exprSet )
QCexpressionMatrix(exprSet, group_list, project = studyID)

library(limma)
design=model.matrix(~ group_list )
fit=lmFit(exprSet,design)
fit=eBayes(fit)
options(digits = 4)

DEG <- topTable(fit,coef=2,adjust='BH',n=Inf)
Volcanic_DEG(DEG)
DEG$symbol <- rownames(DEG)

format_DEG(DEG,studyID)

setwd('test')
input.ds.file = 'input/test.gct'
input.cls.file= 'input/test.cls'
gs.db.file= 'input/c2.MESENCHYMAL.v5.2.symbols.gmt'
output.directory = 'results/'
GSEA(input.ds =  input.ds.file,input.cls =input.cls.file ,gs.db = gs.db.file ,output.directory='./', reshuffling.type      = "gene.labels")


```
