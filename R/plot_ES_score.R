#' re-draw ES-score figures for GSEA
#'
#' It's very complicate,please don't use this unless you have a clear understanding of the mechanism for the GSEA
#' especially for the ES core ,the running ES core, the gene sets !
#' you need to prepare many data by yourself, I really can't help you with them.
#'#'
#' @param Ng  The number of the gene　sets, (eg: 200~300 kegg pathway)
#' @param N   The number of the genes (probably 20,000~30,000 genes)
#' @param phen1  The name for the first phenotype (default:control)
#' @param phen2  The name for the second phenotype (default:case)
#' @param Obs.RES A matrix, running ES score for each gene in specific gene Set.
#' @param Obs.indicator A matrix, 0/1 shows whether a gene belong to a gene Set.
#' @param obs.s2n    A vector of sigal to noise value for each gene (sorted, and Z-score)
#' @param size.G     A vector for the size of each geneSet
#' @param gs.names   A vector for the name of each geneSet
#' @param Obs.ES     A vector for the maximal ES score of each gene Set.
#' @param Obs.arg.ES A vector for the order of the gene shows maximal ES score of each gene Set.
#' @return write some files
#' @export
#' @keywords plot_ES_score
#' @examples
#' #'  plot_ES_score(Ng=12,N=34688,phen1='control',phen2='case',Obs.RES,Obs.indicator,obs.s2n,size.G,gs.names,Obs.ES,Obs.arg.ES,Obs.ES.index)
#'
#'

plot_ES_score <- function(Ng=12,N=34688,phen1='control',phen2='case',Obs.RES,Obs.indicator,obs.s2n,size.G,gs.names,Obs.ES,Obs.arg.ES,Obs.ES.index){
  if(F){
    ## it's commant here,
    setwd('data')
    Obs.RES=read.table('Obs.RES.txt')
    Obs.RES=t(Obs.RES)
    Obs.indicator=read.table('Obs.indicator.txt')
    Obs.indicator=t(Obs.indicator)
    obs.s2n=read.table('obs.s2n.txt')[,1]
    size.G=read.table('size.G.txt')[,1]
    gs.names=read.table('gs.names.txt')[,1]
    Obs.arg.ES=read.table('Obs.arg.ES.txt')[,1]
    Obs.ES.index=read.table('Obs.ES.index.txt')[,1]
    Obs.ES=read.table('Obs.ES.txt')[1,]

  }
  for (i in 1:Ng) {
    png(paste0('geneset_',gs.names[i],'.png'))
    ind <- 1:N
    min.RES <- min(Obs.RES[i,])
    max.RES <- max(Obs.RES[i,])
    if (max.RES < 0.3) max.RES <- 0.3
    if (min.RES > -0.3) min.RES <- -0.3
    delta <- (max.RES - min.RES)*0.50
    min.plot <- min.RES - 2*delta
    max.plot <- max.RES
    max.corr <- max(obs.s2n)
    min.corr <- min(obs.s2n)
    Obs.correl.vector.norm <- (obs.s2n - min.corr)/(max.corr - min.corr)*1.25*delta + min.plot
    zero.corr.line <- (- min.corr/(max.corr - min.corr))*1.25*delta + min.plot
    col <- ifelse(Obs.ES[i] > 0, 2, 4)

    # Running enrichment plot

    sub.string <- paste("Number of genes: ", N, " (in list), ", size.G[i], " (in gene set)", sep = "", collapse="")

    main.string <- paste("Gene Set ", i, ":", gs.names[i])

    plot(ind, Obs.RES[i,], main = main.string, sub = sub.string, xlab = "Gene List Index", ylab = "Running Enrichment Score (RES)", xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = col)
    for (j in seq(1, N, 20)) {
      lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
    }
    lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
    lines(c(Obs.arg.ES[i], Obs.arg.ES[i]), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col) # max enrichment vertical line
    for (j in 1:N) {
      if (Obs.indicator[i, j] == 1) {
        lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
      }
    }
    lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
    lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
    temp <- order(abs(obs.s2n), decreasing=T)
    arg.correl <- temp[N]
    lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 3) # zero crossing correlation vertical line

    leg.txt <- paste("\"", phen1, "\" ", sep="", collapse="")
    text(x=1, y=min.plot, adj = c(0, 0), labels=leg.txt, cex = 1.0)

    leg.txt <- paste("\"", phen2, "\" ", sep="", collapse="")
    text(x=N, y=min.plot, adj = c(1, 0), labels=leg.txt, cex = 1.0)

    adjx <- ifelse(Obs.ES[i] > 0, 0, 1)

    leg.txt <- paste("Peak at ", Obs.arg.ES[i], sep="", collapse="")
    text(x=Obs.arg.ES[i], y=min.plot + 1.8*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)

    leg.txt <- paste("Zero crossing at ", arg.correl, sep="", collapse="")
    text(x=arg.correl, y=min.plot + 1.95*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
    dev.off()
  }

}
