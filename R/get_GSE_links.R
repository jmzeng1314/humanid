#' Get the links for GSE studyID
#'
#'
#' @param studyID A standard study ID
#' @param destdir where to store the files download according to current GSE studyID
#' @return write some files
#' @export
#' @keywords get_GSE_links
#' @examples
#' #' get_GSE_links();get_GSE_links(studyID='GSE1009',down=T);get_GSE_links('GSE42872')
#'
#'

get_GSE_links <- function(studyID='GSE1009',down=F,destdir='./'){
  ## studyID
  ## destdir
  ## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1009/matrix/GSE1009_series_matrix.txt.gz
  ## ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1009/suppl/GSE1009_RAW.tar
  ## http://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&mode=csv&series=1009
  supp_link=paste0(
    "ftp://ftp.ncbi.nlm.nih.gov/geo/series/",
    substr(studyID,1,nchar(studyID)-3),"nnn/",
    studyID,"/suppl/",
    studyID,"_RAW.tar"
  )
  meta_link=paste0(
    "http://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&mode=csv&series=",
    substr(studyID,4,nchar(studyID))
  )
  matrix_link=paste0(
    "ftp://ftp.ncbi.nlm.nih.gov/geo/series/",
    substr(studyID,1,nchar(studyID)-3),"nnn/",
    studyID,"/matrix/",
    studyID,"_series_matrix.txt.gz"
  )
  print (paste0('The URL for raw data is : ',supp_link))
  print (paste0('The URL for metadata is : ',meta_link))
  print (paste0('The URL for matrix   is : ',matrix_link))
  if(down){
    if( !file.exists( file.path(destdir,studyID) ) ){
      dir.create( ( file.path(destdir,studyID) ) )
    }
    #download.file(meta_link,destfile = "./tmp.csv")
    a=read.csv(meta_link)
    write.csv(a,file.path(destdir,paste0(studyID,".meta.txt")))
    download.file(supp_link,
                  destfile = file.path(destdir,studyID,paste0(studyID,"_RAW.tar")),
                  method = 'auto'
    )
    download.file(matrix_link,
                  destfile = file.path(destdir,paste0(studyID,"_series_matrix.txt.gz")),
                  method = 'auto'
    )
    #system('tar xvf *.tar')
    setwd(file.path(destdir,studyID))
    for (i in list.files(pattern = "*tar")){
      untar(i,exdir=destdir)
    }
    all_cel_files=list.files(destdir,pattern="^GSM")
    old_dir <- getwd()
    setwd(destdir)
    for (i in all_cel_files){
      if( grepl("[-_]",i) && (grepl(".CEL.gz",i)  || grepl(".cel.gz",i)  )){
        new=paste(sub("[-_].*","",i),".CEL.gz",sep="")
        file.rename(i,new)
        print(paste0(i,"\t to \t",new,"\n"))
      }
    }
    setwd(old_dir)
    #system('mv *.CEL.gz ../')
  }
}

# get_GSE_links(studyID='GSE1009',down=T)

