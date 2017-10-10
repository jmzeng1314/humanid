# The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT This software and its documentation are copyright 2003
# by the Broad Institute/Massachusetts Institute of Technology.  All rights are reserved.  This software is
# supplied without any warranty or guaranteed support whatsoever. Neither the Broad Institute nor MIT can be
# responsible for its use, misuse, or functionality.


# G S E A -- Gene Set Enrichment Analysis

# Auxiliary functions and definitions

GSEA.GeneRanking <- function(A, class.labels, gene.labels, nperm, permutation.type = 0, sigma.correction = "GeneCluster", 
    fraction = 1, replace = F, reverse.sign = F) {
    
    # This function ranks the genes according to the signal to noise ratio for the actual phenotype and also random
    # permutations and bootstrap subsamples of both the observed and random phenotypes. It uses matrix operations to
    # implement the signal to noise calculation in stages and achieves fast execution speed. It supports two types
    # of permutations: random (unbalanced) and balanced.  It also supports subsampling and bootstrap by using
    # masking and multiple-count variables.  When 'fraction' is set to 1 (default) the there is no subsampling or
    # boostrapping and the matrix of observed signal to noise ratios will have the same value for all permutations.
    # This is wasteful but allows to support all the multiple options with the same code. Notice that the second
    # matrix for the null distribution will still have the values for the random permutations (null distribution).
    # This mode (fraction = 1.0) is the defaults, the recommended one and the one used in the examples.  It is also
    # the one that has be tested more thoroughly. The resampling and boostrapping options are intersting to obtain
    # smooth estimates of the observed distribution but its is left for the expert user who may want to perform some
    # sanity checks before trusting the code.  Inputs: A: Matrix of gene expression values (rows are genes, columns
    # are samples) class.labels: Phenotype of class disticntion of interest. A vector of binary labels having first
    # the 1's and then the 0's gene.labels: gene labels. Vector of probe ids or accession numbers for the rows of
    # the expression matrix nperm: Number of random permutations/bootstraps to perform permutation.type: Permutation
    # type: 0 = unbalanced, 1 = balanced. For experts only (default: 0) sigma.correction: Correction to the signal
    # to noise ratio (Default = GeneCluster, a choice to support the way it was handled in a previous package)
    # fraction: Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0) replace:
    # Resampling mode (replacement or not replacement). For experts only (default: F) reverse.sign: Reverse
    # direction of gene list (default = F) Outputs: s2n.matrix: Matrix with random permuted or bootstraps signal to
    # noise ratios (rows are genes, columns are permutations or bootstrap subsamplings obs.s2n.matrix: Matrix with
    # observed signal to noise ratios (rows are genes, columns are boostraps subsamplings. If fraction is set to 1.0
    # then all the columns have the same values order.matrix: Matrix with the orderings that will sort the columns
    # of the obs.s2n.matrix in decreasing s2n order obs.order.matrix: Matrix with the orderings that will sort the
    # columns of the s2n.matrix in decreasing s2n order The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT This
    # software and its documentation are copyright 2003 by the Broad Institute/Massachusetts Institute of
    # Technology.  All rights are reserved.  This software is supplied without any warranty or guaranteed support
    # whatsoever. Neither the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
    
    A <- A + 1e-08
    
    N <- length(A[, 1])
    Ns <- length(A[1, ])
    
    subset.mask <- matrix(0, nrow = Ns, ncol = nperm)
    reshuffled.class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
    reshuffled.class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
    class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
    class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
    
    order.matrix <- matrix(0, nrow = N, ncol = nperm)
    obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
    s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
    obs.s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
    
    obs.gene.labels <- vector(length = N, mode = "character")
    obs.gene.descs <- vector(length = N, mode = "character")
    obs.gene.symbols <- vector(length = N, mode = "character")
    
    M1 <- matrix(0, nrow = N, ncol = nperm)
    M2 <- matrix(0, nrow = N, ncol = nperm)
    S1 <- matrix(0, nrow = N, ncol = nperm)
    S2 <- matrix(0, nrow = N, ncol = nperm)
    
    gc()
    
    C <- split(class.labels, class.labels)
    class1.size <- length(C[[1]])
    class2.size <- length(C[[2]])
    class1.index <- seq(1, class1.size, 1)
    class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)
    
    for (r in 1:nperm) {
        class1.subset <- sample(class1.index, size = ceiling(class1.size * fraction), replace = replace)
        class2.subset <- sample(class2.index, size = ceiling(class2.size * fraction), replace = replace)
        class1.subset.size <- length(class1.subset)
        class2.subset.size <- length(class2.subset)
        subset.class1 <- rep(0, class1.size)
        for (i in 1:class1.size) {
            if (is.element(class1.index[i], class1.subset)) {
                subset.class1[i] <- 1
            }
        }
        subset.class2 <- rep(0, class2.size)
        for (i in 1:class2.size) {
            if (is.element(class2.index[i], class2.subset)) {
                subset.class2[i] <- 1
            }
        }
        subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
        fraction.class1 <- class1.size/Ns
        fraction.class2 <- class2.size/Ns
        
        if (permutation.type == 0) {
            # random (unbalanced) permutation
            full.subset <- c(class1.subset, class2.subset)
            label1.subset <- sample(full.subset, size = Ns * fraction.class1)
            reshuffled.class.labels1[, r] <- rep(0, Ns)
            reshuffled.class.labels2[, r] <- rep(0, Ns)
            class.labels1[, r] <- rep(0, Ns)
            class.labels2[, r] <- rep(0, Ns)
            for (i in 1:Ns) {
                m1 <- sum(!is.na(match(label1.subset, i)))
                m2 <- sum(!is.na(match(full.subset, i)))
                reshuffled.class.labels1[i, r] <- m1
                reshuffled.class.labels2[i, r] <- m2 - m1
                if (i <= class1.size) {
                  class.labels1[i, r] <- m2
                  class.labels2[i, r] <- 0
                } else {
                  class.labels1[i, r] <- 0
                  class.labels2[i, r] <- m2
                }
            }
            
        } else if (permutation.type == 1) {
            # proportional (balanced) permutation
            
            class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size * fraction.class1))
            class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size * fraction.class1))
            reshuffled.class.labels1[, r] <- rep(0, Ns)
            reshuffled.class.labels2[, r] <- rep(0, Ns)
            class.labels1[, r] <- rep(0, Ns)
            class.labels2[, r] <- rep(0, Ns)
            for (i in 1:Ns) {
                if (i <= class1.size) {
                  m1 <- sum(!is.na(match(class1.label1.subset, i)))
                  m2 <- sum(!is.na(match(class1.subset, i)))
                  reshuffled.class.labels1[i, r] <- m1
                  reshuffled.class.labels2[i, r] <- m2 - m1
                  class.labels1[i, r] <- m2
                  class.labels2[i, r] <- 0
                } else {
                  m1 <- sum(!is.na(match(class2.label1.subset, i)))
                  m2 <- sum(!is.na(match(class2.subset, i)))
                  reshuffled.class.labels1[i, r] <- m1
                  reshuffled.class.labels2[i, r] <- m2 - m1
                  class.labels1[i, r] <- 0
                  class.labels2[i, r] <- m2
                }
            }
        }
    }
    
    # compute S2N for the random permutation matrix
    
    P <- reshuffled.class.labels1 * subset.mask
    n1 <- sum(P[, 1])
    M1 <- A %*% P
    M1 <- M1/n1
    gc()
    A2 <- A * A
    S1 <- A2 %*% P
    S1 <- S1/n1 - M1 * M1
    S1 <- sqrt(abs((n1/(n1 - 1)) * S1))
    gc()
    P <- reshuffled.class.labels2 * subset.mask
    n2 <- sum(P[, 1])
    M2 <- A %*% P
    M2 <- M2/n2
    gc()
    A2 <- A * A
    S2 <- A2 %*% P
    S2 <- S2/n2 - M2 * M2
    S2 <- sqrt(abs((n2/(n2 - 1)) * S2))
    rm(P)
    rm(A2)
    gc()
    
    if (sigma.correction == "GeneCluster") {
        # small sigma 'fix' as used in GeneCluster
        S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
        S2 <- ifelse(S2 == 0, 0.2, S2)
        S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
        S1 <- ifelse(S1 == 0, 0.2, S1)
        gc()
    }
    
    M1 <- M1 - M2
    rm(M2)
    gc()
    S1 <- S1 + S2
    rm(S2)
    gc()
    
    s2n.matrix <- M1/S1
    
    if (reverse.sign == T) {
        s2n.matrix <- -s2n.matrix
    }
    gc()
    
    for (r in 1:nperm) {
        order.matrix[, r] <- order(s2n.matrix[, r], decreasing = T)
    }
    
    # compute S2N for the 'observed' permutation matrix
    
    P <- class.labels1 * subset.mask
    n1 <- sum(P[, 1])
    M1 <- A %*% P
    M1 <- M1/n1
    gc()
    A2 <- A * A
    S1 <- A2 %*% P
    S1 <- S1/n1 - M1 * M1
    S1 <- sqrt(abs((n1/(n1 - 1)) * S1))
    gc()
    P <- class.labels2 * subset.mask
    n2 <- sum(P[, 1])
    M2 <- A %*% P
    M2 <- M2/n2
    gc()
    A2 <- A * A
    S2 <- A2 %*% P
    S2 <- S2/n2 - M2 * M2
    S2 <- sqrt(abs((n2/(n2 - 1)) * S2))
    rm(P)
    rm(A2)
    gc()
    
    if (sigma.correction == "GeneCluster") {
        # small sigma 'fix' as used in GeneCluster
        S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
        S2 <- ifelse(S2 == 0, 0.2, S2)
        S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
        S1 <- ifelse(S1 == 0, 0.2, S1)
        gc()
    }
    
    M1 <- M1 - M2
    rm(M2)
    gc()
    S1 <- S1 + S2
    rm(S2)
    gc()
    
    obs.s2n.matrix <- M1/S1
    gc()
    
    if (reverse.sign == T) {
        obs.s2n.matrix <- -obs.s2n.matrix
    }
    
    for (r in 1:nperm) {
        obs.order.matrix[, r] <- order(obs.s2n.matrix[, r], decreasing = T)
    }
    
    return(list(s2n.matrix = s2n.matrix, obs.s2n.matrix = obs.s2n.matrix, order.matrix = order.matrix, obs.order.matrix = obs.order.matrix))
}

GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {
    # Computes the weighted GSEA score of gene.set in gene.list.  The weighted score type is the exponent of the
    # correlation weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score
    # type is 1 or 2 it is necessary to input the correlation vector with the values in the same order as in the
    # gene list.  Inputs: gene.list: The ordered gene list (e.g. integers indicating the original position in the
    # input dataset) gene.set: A gene set (e.g. integers indicating the location of those genes in the input
    # dataset) weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2
    # (over-weighted) correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to
    # the genes in the gene list Outputs: ES: Enrichment score (real number between -1 and +1) arg.ES: Location in
    # gene.list where the peak running enrichment occurs (peak of the 'mountain') RES: Numerical vector containing
    # the running enrichment score for all locations in the gene list tag.indicator: Binary vector indicating the
    # location of the gene sets (1's) in the gene list The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT This
    # software and its documentation are copyright 2003 by the Broad Institute/Massachusetts Institute of
    # Technology.  All rights are reserved.  This software is supplied without any warranty or guaranteed support
    # whatsoever. Neither the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
    
    tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))  # notice that the sign is 0 (no tag) or 1 (tag)
    no.tag.indicator <- 1 - tag.indicator
    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <- N - Nh
    if (weighted.score.type == 0) {
        correl.vector <- rep(1, N)
    }
    alpha <- weighted.score.type
    correl.vector <- abs(correl.vector^alpha)
    sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
    norm.tag <- 1/sum.correl.tag
    norm.no.tag <- 1/Nm
    RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
        # ES <- max.ES
        ES <- signif(max.ES, digits = 5)
        arg.ES <- which.max(RES)
    } else {
        # ES <- min.ES
        ES <- signif(min.ES, digits = 5)
        arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}


OLD.GSEA.EnrichmentScore <- function(gene.list, gene.set) {
    # Computes the original GSEA score from Mootha et al 2003 of gene.set in gene.list Inputs: gene.list: The
    # ordered gene list (e.g. integers indicating the original position in the input dataset) gene.set: A gene set
    # (e.g. integers indicating the location of those genes in the input dataset) Outputs: ES: Enrichment score
    # (real number between -1 and +1) arg.ES: Location in gene.list where the peak running enrichment occurs (peak
    # of the 'mountain') RES: Numerical vector containing the running enrichment score for all locations in the gene
    # list tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list The Broad
    # Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT This software and its documentation are copyright 2003 by the
    # Broad Institute/Massachusetts Institute of Technology.  All rights are reserved.  This software is supplied
    # without any warranty or guaranteed support whatsoever. Neither the Broad Institute nor MIT can be responsible
    # for its use, misuse, or functionality.
    
    tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))  # notice that the sign is 0 (no tag) or 1 (tag)
    no.tag.indicator <- 1 - tag.indicator
    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <- N - Nh
    
    norm.tag <- sqrt((N - Nh)/Nh)
    norm.no.tag <- sqrt(Nh/(N - Nh))
    
    RES <- cumsum(tag.indicator * norm.tag - no.tag.indicator * norm.no.tag)
    max.ES <- max(RES)
    min.ES <- min(RES)
    if (max.ES > -min.ES) {
        ES <- signif(max.ES, digits = 5)
        arg.ES <- which.max(RES)
    } else {
        ES <- signif(min.ES, digits = 5)
        arg.ES <- which.min(RES)
    }
    return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}

GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {
    # Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in
    # GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.  This call
    # is intended to be used to asses the enrichment of random permutations rather than the observed one.  The
    # weighted score type is the exponent of the correlation weight: 0 (unweighted = Kolmogorov-Smirnov), 1
    # (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is necessary to input the correlation
    # vector with the values in the same order as in the gene list.  Inputs: gene.list: The ordered gene list (e.g.
    # integers indicating the original position in the input dataset) gene.set: A gene set (e.g. integers indicating
    # the location of those genes in the input dataset) weighted.score.type: Type of score: weight: 0 (unweighted =
    # Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted) correl.vector: A vector with the coorelations (e.g.
    # signal to noise scores) corresponding to the genes in the gene list Outputs: ES: Enrichment score (real number
    # between -1 and +1) The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT This software and its documentation
    # are copyright 2003 by the Broad Institute/Massachusetts Institute of Technology.  All rights are reserved.
    # This software is supplied without any warranty or guaranteed support whatsoever. Neither the Broad Institute
    # nor MIT can be responsible for its use, misuse, or functionality.
    
    N <- length(gene.list)
    Nh <- length(gene.set)
    Nm <- N - Nh
    
    loc.vector <- vector(length = N, mode = "numeric")
    peak.res.vector <- vector(length = Nh, mode = "numeric")
    valley.res.vector <- vector(length = Nh, mode = "numeric")
    tag.correl.vector <- vector(length = Nh, mode = "numeric")
    tag.diff.vector <- vector(length = Nh, mode = "numeric")
    tag.loc.vector <- vector(length = Nh, mode = "numeric")
    
    loc.vector[gene.list] <- seq(1, N)
    tag.loc.vector <- loc.vector[gene.set]
    
    tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
    
    if (weighted.score.type == 0) {
        tag.correl.vector <- rep(1, Nh)
    } else if (weighted.score.type == 1) {
        tag.correl.vector <- correl.vector[tag.loc.vector]
        tag.correl.vector <- abs(tag.correl.vector)
    } else if (weighted.score.type == 2) {
        tag.correl.vector <- correl.vector[tag.loc.vector] * correl.vector[tag.loc.vector]
        tag.correl.vector <- abs(tag.correl.vector)
    } else {
        tag.correl.vector <- correl.vector[tag.loc.vector]^weighted.score.type
        tag.correl.vector <- abs(tag.correl.vector)
    }
    
    norm.tag <- 1/sum(tag.correl.vector)
    tag.correl.vector <- tag.correl.vector * norm.tag
    norm.no.tag <- 1/Nm
    tag.diff.vector[1] <- (tag.loc.vector[1] - 1)
    tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
    tag.diff.vector <- tag.diff.vector * norm.no.tag
    peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
    valley.res.vector <- peak.res.vector - tag.correl.vector
    max.ES <- max(peak.res.vector)
    min.ES <- min(valley.res.vector)
    ES <- signif(ifelse(max.ES > -min.ES, max.ES, min.ES), digits = 5)
    
    return(list(ES = ES))
    
}

GSEA.HeatMapPlot <- function(V, row.names = F, col.labels, col.classes, col.names = F, main = " ", xlab = " ", ylab = " ") {
    # Plots a heatmap 'pinkogram' of a gene expression matrix including phenotype vector and gene, sample and
    # phenotype labels The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT This software and its documentation
    # are copyright 2003 by the Broad Institute/Massachusetts Institute of Technology.  All rights are reserved.
    # This software is supplied without any warranty or guaranteed support whatsoever. Neither the Broad Institute
    # nor MIT can be responsible for its use, misuse, or functionality.
    
    n.rows <- length(V[, 1])
    n.cols <- length(V[1, ])
    row.mean <- apply(V, MARGIN = 1, FUN = mean)
    row.sd <- apply(V, MARGIN = 1, FUN = sd)
    row.n <- length(V[, 1])
    for (i in 1:n.rows) {
        if (row.sd[i] == 0) {
            V[i, ] <- 0
        } else {
            V[i, ] <- (V[i, ] - row.mean[i])/(0.5 * row.sd[i])
        }
        V[i, ] <- ifelse(V[i, ] < -6, -6, V[i, ])
        V[i, ] <- ifelse(V[i, ] > 6, 6, V[i, ])
    }
    
    mycol <- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", 
        "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000")  # blue-pinkogram colors. The first and last are the colors to indicate the class vector (phenotype). This is the 1998-vintage, pre-gene cluster, original pinkogram color map
    
    mid.range.V <- mean(range(V)) - 0.1
    heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
    heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]
    heatm[n.rows + 1, ] <- ifelse(col.labels == 0, 7, -7)
    image(1:n.cols, 1:(n.rows + 1), t(heatm), col = mycol, axes = FALSE, main = main, xlab = xlab, ylab = ylab)
    
    if (length(row.names) > 1) {
        numC <- nchar(row.names)
        size.row.char <- 35/(n.rows + 5)
        size.col.char <- 25/(n.cols + 5)
        maxl <- floor(n.rows/1.6)
        for (i in 1:n.rows) {
            row.names[i] <- substr(row.names[i], 1, maxl)
        }
        row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
        axis(2, at = 1:(n.rows + 1), labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char, 
            font.axis = 2, line = -1)
    }
    
    if (length(col.names) > 1) {
        axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, 
            line = -1)
    }
    
    C <- split(col.labels, col.labels)
    class1.size <- length(C[[1]])
    class2.size <- length(C[[2]])
    axis(3, at = c(floor(class1.size/2), class1.size + floor(class2.size/2)), labels = col.classes, tick = FALSE, 
        las = 1, cex.axis = 1.25, font.axis = 2, line = -1)
    
    return()
}

GSEA.Res2Frame <- function(filename = "NULL") {
    # Reads a gene expression dataset in RES format and converts it into an R data frame The Broad Institute
    # SOFTWARE COPYRIGHT NOTICE AGREEMENT This software and its documentation are copyright 2003 by the Broad
    # Institute/Massachusetts Institute of Technology.  All rights are reserved.  This software is supplied without
    # any warranty or guaranteed support whatsoever. Neither the Broad Institute nor MIT can be responsible for its
    # use, misuse, or functionality.
    
    header.cont <- readLines(filename, n = 1)
    temp <- unlist(strsplit(header.cont, "\t"))
    colst <- length(temp)
    header.labels <- temp[seq(3, colst, 2)]
    ds <- read.delim(filename, header = F, row.names = 2, sep = "\t", skip = 3, blank.lines.skip = T, comment.char = "", 
        as.is = T)
    colst <- length(ds[1, ])
    cols <- (colst - 1)/2
    rows <- length(ds[, 1])
    A <- matrix(nrow = rows - 1, ncol = cols)
    A <- ds[1:rows, seq(2, colst, 2)]
    table1 <- data.frame(A)
    names(table1) <- header.labels
    return(table1)
}

GSEA.Gct2Frame <- function(filename = "NULL") {
    # Reads a gene expression dataset in GCT format and converts it into an R data frame The Broad Institute
    # SOFTWARE COPYRIGHT NOTICE AGREEMENT This software and its documentation are copyright 2003 by the Broad
    # Institute/Massachusetts Institute of Technology.  All rights are reserved.  This software is supplied without
    # any warranty or guaranteed support whatsoever. Neither the Broad Institute nor MIT can be responsible for its
    # use, misuse, or functionality.
    ds <- read.delim(filename, header = T, sep = "\t", skip = 2, row.names = 1, blank.lines.skip = T, comment.char = "", 
        as.is = T)
    ds <- ds[-1]
    return(ds)
}

GSEA.Gct2Frame2 <- function(filename = "NULL") {
    # Reads a gene expression dataset in GCT format and converts it into an R data frame The Broad Institute
    # SOFTWARE COPYRIGHT NOTICE AGREEMENT This software and its documentation are copyright 2003 by the Broad
    # Institute/Massachusetts Institute of Technology.  All rights are reserved.  This software is supplied without
    # any warranty or guaranteed support whatsoever. Neither the Broad Institute nor MIT can be responsible for its
    # use, misuse, or functionality.
    content <- readLines(filename)
    content <- content[-1]
    content <- content[-1]
    col.names <- noquote(unlist(strsplit(content[1], "\t")))
    col.names <- col.names[c(-1, -2)]
    num.cols <- length(col.names)
    content <- content[-1]
    num.lines <- length(content)
    
    
    row.nam <- vector(length = num.lines, mode = "character")
    row.des <- vector(length = num.lines, mode = "character")
    m <- matrix(0, nrow = num.lines, ncol = num.cols)
    
    for (i in 1:num.lines) {
        line.list <- noquote(unlist(strsplit(content[i], "\t")))
        row.nam[i] <- noquote(line.list[1])
        row.des[i] <- noquote(line.list[2])
        line.list <- line.list[c(-1, -2)]
        for (j in 1:length(line.list)) {
            m[i, j] <- as.numeric(line.list[j])
        }
    }
    ds <- data.frame(m)
    names(ds) <- col.names
    row.names(ds) <- row.nam
    return(ds)
}

GSEA.ReadClsFile <- function(file = "NULL") {
    # Reads a class vector CLS file and defines phenotype and class labels vectors for the samples in a gene
    # expression file (RES or GCT format) The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT This software and
    # its documentation are copyright 2003 by the Broad Institute/Massachusetts Institute of Technology.  All rights
    # are reserved.  This software is supplied without any warranty or guaranteed support whatsoever. Neither the
    # Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
    
    cls.cont <- readLines(file)
    num.lines <- length(cls.cont)
    class.list <- unlist(strsplit(cls.cont[[3]], " "))
    s <- length(class.list)
    t <- table(class.list)
    l <- length(t)
    phen <- vector(length = l, mode = "character")
    phen.label <- vector(length = l, mode = "numeric")
    class.v <- vector(length = s, mode = "numeric")
    for (i in 1:l) {
        phen[i] <- noquote(names(t)[i])
        phen.label[i] <- i - 1
    }
    for (i in 1:s) {
        for (j in 1:l) {
            if (class.list[i] == phen[j]) {
                class.v[i] <- phen.label[j]
            }
        }
    }
    return(list(phen = phen, class.v = class.v))
}

GSEA.Threshold <- function(V, thres, ceil) {
    # Threshold and ceiling pre-processing for gene expression matrix The Broad Institute SOFTWARE COPYRIGHT NOTICE
    # AGREEMENT This software and its documentation are copyright 2003 by the Broad Institute/Massachusetts
    # Institute of Technology.  All rights are reserved.  This software is supplied without any warranty or
    # guaranteed support whatsoever. Neither the Broad Institute nor MIT can be responsible for its use, misuse, or
    # functionality.
    
    V[V < thres] <- thres
    V[V > ceil] <- ceil
    return(V)
}

GSEA.VarFilter <- function(V, fold, delta, gene.names = "NULL") {
    # Variation filter pre-processing for gene expression matrix The Broad Institute SOFTWARE COPYRIGHT NOTICE
    # AGREEMENT This software and its documentation are copyright 2003 by the Broad Institute/Massachusetts
    # Institute of Technology.  All rights are reserved.  This software is supplied without any warranty or
    # guaranteed support whatsoever. Neither the Broad Institute nor MIT can be responsible for its use, misuse, or
    # functionality.
    
    cols <- length(V[1, ])
    rows <- length(V[, 1])
    row.max <- apply(V, MARGIN = 1, FUN = max)
    row.min <- apply(V, MARGIN = 1, FUN = min)
    flag <- array(dim = rows)
    flag <- (row.max/row.min > fold) & (row.max - row.min > delta)
    size <- sum(flag)
    B <- matrix(0, nrow = size, ncol = cols)
    j <- 1
    if (gene.names == "NULL") {
        for (i in 1:rows) {
            if (flag[i]) {
                B[j, ] <- V[i, ]
                j <- j + 1
            }
        }
        return(B)
    } else {
        new.list <- vector(mode = "character", length = size)
        for (i in 1:rows) {
            if (flag[i]) {
                B[j, ] <- V[i, ]
                new.list[j] <- gene.names[i]
                j <- j + 1
            }
        }
        return(list(V = B, new.list = new.list))
    }
}

GSEA.NormalizeRows <- function(V) {
    # Stardardize rows of a gene expression matrix The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT This
    # software and its documentation are copyright 2003 by the Broad Institute/Massachusetts Institute of
    # Technology.  All rights are reserved.  This software is supplied without any warranty or guaranteed support
    # whatsoever. Neither the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
    
    row.mean <- apply(V, MARGIN = 1, FUN = mean)
    row.sd <- apply(V, MARGIN = 1, FUN = sd)
    row.n <- length(V[, 1])
    for (i in 1:row.n) {
        if (row.sd[i] == 0) {
            V[i, ] <- 0
        } else {
            V[i, ] <- (V[i, ] - row.mean[i])/row.sd[i]
        }
    }
    return(V)
}

GSEA.NormalizeCols <- function(V) {
    # Stardardize columns of a gene expression matrix The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT This
    # software and its documentation are copyright 2003 by the Broad Institute/Massachusetts Institute of
    # Technology.  All rights are reserved.  This software is supplied without any warranty or guaranteed support
    # whatsoever. Neither the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
    
    col.mean <- apply(V, MARGIN = 2, FUN = mean)
    col.sd <- apply(V, MARGIN = 2, FUN = sd)
    col.n <- length(V[1, ])
    for (i in 1:col.n) {
        if (col.sd[i] == 0) {
            V[i, ] <- 0
        } else {
            V[, i] <- (V[, i] - col.mean[i])/col.sd[i]
        }
    }
    return(V)
}

# end of auxiliary functions

# ----------------------------------------------------------------------------------------


GSEA.ConsPlot <- function(V, col.names, main = " ", sub = " ", xlab = " ", ylab = " ") {
    
    # Plots a heatmap plot of a consensus matrix
    
    cols <- length(V[1, ])
    B <- matrix(0, nrow = cols, ncol = cols)
    max.val <- max(V)
    min.val <- min(V)
    for (i in 1:cols) {
        for (j in 1:cols) {
            k <- cols - i + 1
            B[k, j] <- max.val - V[i, j] + min.val
        }
    }
    
    
    
    # col.map <- c(rainbow(100, s = 1.0, v = 0.75, start = 0.0, end = 0.75, gamma = 1.5), '#BBBBBB', '#333333',
    # '#FFFFFF')
    col.map <- rev(c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", 
        "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D"))
    
    # max.size <- max(nchar(col.names))
    par(mar = c(5, 15, 15, 5))
    image(1:cols, 1:cols, t(B), col = col.map, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)
    
    for (i in 1:cols) {
        col.names[i] <- substr(col.names[i], 1, 25)
    }
    col.names2 <- rev(col.names)
    
    size.col.char <- ifelse(cols < 15, 1, sqrt(15/cols))
    
    axis(2, at = 1:cols, labels = col.names2, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.col.char, font.axis = 1, 
        line = -1)
    axis(3, at = 1:cols, labels = col.names, adj = 1, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 1, 
        line = -1)
    
    return()
}

GSEA.HeatMapPlot2 <- function(V, row.names = "NA", col.names = "NA", main = " ", sub = " ", xlab = " ", ylab = " ", 
    color.map = "default") {
    # Plots a heatmap of a matrix
    
    n.rows <- length(V[, 1])
    n.cols <- length(V[1, ])
    
    if (color.map == "default") {
        color.map <- rev(rainbow(100, s = 1, v = 0.75, start = 0, end = 0.75, gamma = 1.5))
    }
    
    heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
    heatm[1:n.rows, ] <- V[seq(n.rows, 1, -1), ]
    
    par(mar = c(7, 15, 5, 5))
    image(1:n.cols, 1:n.rows, t(heatm), col = color.map, axes = FALSE, main = main, sub = sub, xlab = xlab, ylab = ylab)
    
    if (length(row.names) > 1) {
        size.row.char <- ifelse(n.rows < 15, 1, sqrt(15/n.rows))
        size.col.char <- ifelse(n.cols < 15, 1, sqrt(10/n.cols))
        # size.col.char <- ifelse(n.cols < 2.5, 1, sqrt(2.5/n.cols))
        for (i in 1:n.rows) {
            row.names[i] <- substr(row.names[i], 1, 40)
        }
        row.names <- row.names[seq(n.rows, 1, -1)]
        axis(2, at = 1:n.rows, labels = row.names, adj = 0.5, tick = FALSE, las = 1, cex.axis = size.row.char, font.axis = 1, 
            line = -1)
    }
    
    if (length(col.names) > 1) {
        axis(1, at = 1:n.cols, labels = col.names, tick = FALSE, las = 3, cex.axis = size.col.char, font.axis = 2, 
            line = -1)
    }
    
    return()
}


GSEA.Analyze.Sets <- function(directory, topgs = "", non.interactive.run = F, height = 12, width = 17) {
    
    file.list <- list.files(directory)
    files <- file.list[regexpr(pattern = ".report.", file.list) > 1]
    max.sets <- length(files)
    
    set.table <- matrix(nrow = max.sets, ncol = 5)
    
    for (i in 1:max.sets) {
        temp1 <- strsplit(files[i], split = ".report.")
        temp2 <- strsplit(temp1[[1]][1], split = ".")
        s <- length(temp2[[1]])
        prefix.name <- paste(temp2[[1]][1:(s - 1)], sep = "", collapse = "")
        set.name <- temp2[[1]][s]
        temp3 <- strsplit(temp1[[1]][2], split = ".")
        phenotype <- temp3[[1]][1]
        seq.number <- temp3[[1]][2]
        dataset <- paste(temp2[[1]][1:(s - 1)], sep = "", collapse = ".")
        
        set.table[i, 1] <- files[i]
        
        set.table[i, 3] <- phenotype
        set.table[i, 4] <- as.numeric(seq.number)
        set.table[i, 5] <- dataset
        
        # set.table[i, 2] <- paste(set.name, dataset, sep ='', collapse='')
        set.table[i, 2] <- substr(set.name, 1, 20)
    }
    
    print(c("set name=", prefix.name))
    doc.string <- prefix.name
    
    set.table <- noquote(set.table)
    phen.order <- order(set.table[, 3], decreasing = T)
    set.table <- set.table[phen.order, ]
    phen1 <- names(table(set.table[, 3]))[1]
    phen2 <- names(table(set.table[, 3]))[2]
    set.table.phen1 <- set.table[set.table[, 3] == phen1, ]
    set.table.phen2 <- set.table[set.table[, 3] == phen2, ]
    
    seq.order <- order(as.numeric(set.table.phen1[, 4]), decreasing = F)
    set.table.phen1 <- set.table.phen1[seq.order, ]
    seq.order <- order(as.numeric(set.table.phen2[, 4]), decreasing = F)
    set.table.phen2 <- set.table.phen2[seq.order, ]
    
    # max.sets.phen1 <- length(set.table.phen1[,1]) max.sets.phen2 <- length(set.table.phen2[,1])
    
    if (topgs == "") {
        max.sets.phen1 <- length(set.table.phen1[, 1])
        max.sets.phen2 <- length(set.table.phen2[, 1])
    } else {
        max.sets.phen1 <- ifelse(topgs > length(set.table.phen1[, 1]), length(set.table.phen1[, 1]), topgs)
        max.sets.phen2 <- ifelse(topgs > length(set.table.phen2[, 1]), length(set.table.phen2[, 1]), topgs)
    }
    
    # Analysis for phen1
    
    leading.lists <- NULL
    for (i in 1:max.sets.phen1) {
        inputfile <- paste(directory, set.table.phen1[i, 1], sep = "", collapse = "")
        gene.set <- read.table(file = inputfile, sep = "\t", header = T, comment.char = "", as.is = T)
        leading.set <- as.vector(gene.set[gene.set[, "CORE_ENRICHMENT"] == "YES", "SYMBOL"])
        leading.lists <- c(leading.lists, list(leading.set))
        if (i == 1) {
            all.leading.genes <- leading.set
        } else {
            all.leading.genes <- union(all.leading.genes, leading.set)
        }
    }
    max.genes <- length(all.leading.genes)
    M <- matrix(0, nrow = max.sets.phen1, ncol = max.genes)
    for (i in 1:max.sets.phen1) {
        M[i, ] <- sign(match(all.leading.genes, as.vector(leading.lists[[i]]), nomatch = 0))  # notice that the sign is 0 (no tag) or 1 (tag)
    }
    
    Inter <- matrix(0, nrow = max.sets.phen1, ncol = max.sets.phen1)
    for (i in 1:max.sets.phen1) {
        for (j in 1:max.sets.phen1) {
            Inter[i, j] <- length(intersect(leading.lists[[i]], leading.lists[[j]]))/length(union(leading.lists[[i]], 
                leading.lists[[j]]))
        }
    }
    
    Itable <- data.frame(Inter)
    names(Itable) <- set.table.phen1[1:max.sets.phen1, 2]
    row.names(Itable) <- set.table.phen1[1:max.sets.phen1, 2]
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.overlap.", phen1, sep = "", collapse = "")
            windows(height = width, width = width)
        } else if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.overlap.", phen1, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = width, width = width)
        }
    } else {
        if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.overlap.", phen1, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = width, width = width)
        } else if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.overlap.", phen1, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = width, width = width)
        }
    }
    
    GSEA.ConsPlot(Itable, col.names = set.table.phen1[1:max.sets.phen1, 2], main = " ", sub = paste("Leading Subsets Overlap ", 
        doc.string, " - ", phen1, sep = ""), xlab = " ", ylab = " ")
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = filename, type = "jpeg", device = dev.cur())
            dev.off()
        } else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    } else {
        dev.off()
    }
    
    # Save leading subsets in a GCT file
    
    D.phen1 <- data.frame(M)
    names(D.phen1) <- all.leading.genes
    row.names(D.phen1) <- set.table.phen1[1:max.sets.phen1, 2]
    output <- paste(directory, doc.string, ".leading.genes.", phen1, ".gct", sep = "")
    GSEA.write.gct(D.phen1, filename = output)
    
    # Save leading subsets as a single gene set in a .gmt file
    
    row.header <- paste(doc.string, ".all.leading.genes.", phen1, sep = "")
    output.line <- paste(all.leading.genes, sep = "\t", collapse = "\t")
    output.line <- paste(row.header, row.header, output.line, sep = "\t", collapse = "")
    output <- paste(directory, doc.string, ".all.leading.genes.", phen1, ".gmt", sep = "")
    write(noquote(output.line), file = output, ncolumns = length(output.line))
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.assignment.", phen1, sep = "", collapse = "")
            windows(height = height, width = width)
        } else if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.assignment.", phen1, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = height, width = width)
        }
    } else {
        if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.assignment.", phen1, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = height, width = width)
        } else if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.assignment.", phen1, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = height, width = width)
        }
    }
    
    cmap <- c("#AAAAFF", "#111166")
    GSEA.HeatMapPlot2(V = data.matrix(D.phen1), row.names = row.names(D.phen1), col.names = names(D.phen1), main = "Leading Subsets Assignment", 
        sub = paste(doc.string, " - ", phen1, sep = ""), xlab = " ", ylab = " ", color.map = cmap)
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = filename, type = "jpeg", device = dev.cur())
            dev.off()
        } else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    } else {
        dev.off()
    }
    
    DT1.phen1 <- data.matrix(t(D.phen1))
    DT2.phen1 <- data.frame(DT1.phen1)
    names(DT2.phen1) <- set.table.phen1[1:max.sets.phen1, 2]
    row.names(DT2.phen1) <- all.leading.genes
    # GSEA.write.gct(DT2.phen1, filename=outputfile2.phen1)
    
    # Analysis for phen2
    
    leading.lists <- NULL
    for (i in 1:max.sets.phen2) {
        inputfile <- paste(directory, set.table.phen2[i, 1], sep = "", collapse = "")
        gene.set <- read.table(file = inputfile, sep = "\t", header = T, comment.char = "", as.is = T)
        leading.set <- as.vector(gene.set[gene.set[, "CORE_ENRICHMENT"] == "YES", "SYMBOL"])
        leading.lists <- c(leading.lists, list(leading.set))
        if (i == 1) {
            all.leading.genes <- leading.set
        } else {
            all.leading.genes <- union(all.leading.genes, leading.set)
        }
    }
    max.genes <- length(all.leading.genes)
    M <- matrix(0, nrow = max.sets.phen2, ncol = max.genes)
    for (i in 1:max.sets.phen2) {
        M[i, ] <- sign(match(all.leading.genes, as.vector(leading.lists[[i]]), nomatch = 0))  # notice that the sign is 0 (no tag) or 1 (tag)
    }
    
    Inter <- matrix(0, nrow = max.sets.phen2, ncol = max.sets.phen2)
    for (i in 1:max.sets.phen2) {
        for (j in 1:max.sets.phen2) {
            Inter[i, j] <- length(intersect(leading.lists[[i]], leading.lists[[j]]))/length(union(leading.lists[[i]], 
                leading.lists[[j]]))
        }
    }
    
    Itable <- data.frame(Inter)
    names(Itable) <- set.table.phen2[1:max.sets.phen2, 2]
    row.names(Itable) <- set.table.phen2[1:max.sets.phen2, 2]
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.overlap.", phen2, sep = "", collapse = "")
            windows(height = width, width = width)
        } else if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.overlap.", phen2, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = width, width = width)
        }
    } else {
        if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.overlap.", phen2, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = width, width = width)
        } else if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.overlap.", phen2, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = width, width = width)
        }
    }
    
    GSEA.ConsPlot(Itable, col.names = set.table.phen2[1:max.sets.phen2, 2], main = " ", sub = paste("Leading Subsets Overlap ", 
        doc.string, " - ", phen2, sep = ""), xlab = " ", ylab = " ")
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = filename, type = "jpeg", device = dev.cur())
            dev.off()
        } else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    } else {
        dev.off()
    }
    
    # Save leading subsets in a GCT file
    
    D.phen2 <- data.frame(M)
    names(D.phen2) <- all.leading.genes
    row.names(D.phen2) <- set.table.phen2[1:max.sets.phen2, 2]
    output <- paste(directory, doc.string, ".leading.genes.", phen2, ".gct", sep = "")
    GSEA.write.gct(D.phen2, filename = output)
    
    # Save primary subsets as a single gene set in a .gmt file
    
    row.header <- paste(doc.string, ".all.leading.genes.", phen2, sep = "")
    output.line <- paste(all.leading.genes, sep = "\t", collapse = "\t")
    output.line <- paste(row.header, row.header, output.line, sep = "\t", collapse = "")
    output <- paste(directory, doc.string, ".all.leading.genes.", phen2, ".gmt", sep = "")
    write(noquote(output.line), file = output, ncolumns = length(output.line))
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.assignment.", phen2, sep = "", collapse = "")
            windows(height = height, width = width)
        } else if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.assignment.", phen2, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = height, width = width)
        }
    } else {
        if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.assignment.", phen2, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = height, width = width)
        } else if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.assignment.", phen2, ".pdf", sep = "", collapse = "")
            pdf(file = filename, height = height, width = width)
        }
    }
    
    cmap <- c("#AAAAFF", "#111166")
    GSEA.HeatMapPlot2(V = data.matrix(D.phen2), row.names = row.names(D.phen2), col.names = names(D.phen2), main = "Leading Subsets Assignment", 
        sub = paste(doc.string, " - ", phen2, sep = ""), xlab = " ", ylab = " ", color.map = cmap)
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = filename, type = "jpeg", device = dev.cur())
            dev.off()
        } else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    } else {
        dev.off()
    }
    
    DT1.phen2 <- data.matrix(t(D.phen2))
    DT2.phen2 <- data.frame(DT1.phen2)
    names(DT2.phen2) <- set.table.phen2[1:max.sets.phen2, 2]
    row.names(DT2.phen2) <- all.leading.genes
    # GSEA.write.gct(DT2.phen2, filename=outputfile2.phen2)
    
    # Resort columns and rows for phen1
    
    A <- data.matrix(D.phen1)
    A.row.names <- row.names(D.phen1)
    A.names <- names(D.phen1)
    
    # Max.genes
    
    # init <- 1 for (k in 1:max.sets.phen1) { end <- which.max(cumsum(A[k,])) if (end - init > 1) { B <-
    # A[,init:end] B.names <- A.names[init:end] dist.matrix <- dist(t(B)) HC <- hclust(dist.matrix,
    # method='average') B <- B[,HC$order] + 0.2*(k %% 2) B <- B[,HC$order] A[,init:end] <- B A.names[init:end] <-
    # B.names[HC$order] init <- end + 1 } }
    
    # windows(width=14, height=10) GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, sub = ' ',
    # main = paste('Primary Sets Assignment - ', doc.string, ' - ', phen1, sep=''), xlab=' ', ylab=' ')
    
    dist.matrix <- dist(t(A))
    HC <- hclust(dist.matrix, method = "average")
    A <- A[, HC$order]
    A.names <- A.names[HC$order]
    
    dist.matrix <- dist(A)
    HC <- hclust(dist.matrix, method = "average")
    A <- A[HC$order, ]
    A.row.names <- A.row.names[HC$order]
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, sep = "", collapse = "")
            windows(height = height, width = width)
        } else if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep = "", 
                collapse = "")
            pdf(file = filename, height = height, width = width)
        }
    } else {
        if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep = "", 
                collapse = "")
            pdf(file = filename, height = height, width = width)
        } else if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, ".pdf", sep = "", 
                collapse = "")
            pdf(file = filename, height = height, width = width)
        }
    }
    
    
    
    cmap <- c("#AAAAFF", "#111166")
    # GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, main = 'Leading Subsets Assignment
    # (clustered)', sub = paste(doc.string, ' - ', phen1, sep=''), xlab=' ', ylab=' ', color.map = cmap)
    
    GSEA.HeatMapPlot2(V = t(A), row.names = A.names, col.names = A.row.names, main = "Leading Subsets Assignment (clustered)", 
        sub = paste(doc.string, " - ", phen1, sep = ""), xlab = " ", ylab = " ", color.map = cmap)
    
    text.filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen1, ".txt", sep = "", collapse = "")
    line.list <- c("Gene", A.row.names)
    line.header <- paste(line.list, collapse = "\t")
    line.length <- length(A.row.names) + 1
    write(line.header, file = text.filename, ncolumns = line.length)
    write.table(t(A), file = text.filename, append = T, quote = F, col.names = F, row.names = T, sep = "\t")
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = filename, type = "jpeg", device = dev.cur())
            dev.off()
        } else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    } else {
        dev.off()
    }
    
    
    
    
    
    
    # resort columns and rows for phen2
    
    A <- data.matrix(D.phen2)
    A.row.names <- row.names(D.phen2)
    A.names <- names(D.phen2)
    
    # Max.genes
    
    # init <- 1 for (k in 1:max.sets.phen2) { end <- which.max(cumsum(A[k,])) if (end - init > 1) { B <-
    # A[,init:end] B.names <- A.names[init:end] dist.matrix <- dist(t(B)) HC <- hclust(dist.matrix,
    # method='average') B <- B[,HC$order] + 0.2*(k %% 2) B <- B[,HC$order] A[,init:end] <- B A.names[init:end] <-
    # B.names[HC$order] init <- end + 1 } }
    
    # windows(width=14, height=10) GESA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, sub = ' ',
    # main = paste('Primary Sets Assignment - ', doc.string, ' - ', phen2, sep=''), xlab=' ', ylab=' ')
    
    dist.matrix <- dist(t(A))
    HC <- hclust(dist.matrix, method = "average")
    A <- A[, HC$order]
    A.names <- A.names[HC$order]
    
    dist.matrix <- dist(A)
    HC <- hclust(dist.matrix, method = "average")
    A <- A[HC$order, ]
    A.row.names <- A.row.names[HC$order]
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, sep = "", collapse = "")
            windows(height = height, width = width)
        } else if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep = "", 
                collapse = "")
            pdf(file = filename, height = height, width = width)
        }
    } else {
        if (.Platform$OS.type == "unix") {
            filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep = "", 
                collapse = "")
            pdf(file = filename, height = height, width = width)
        } else if (.Platform$OS.type == "windows") {
            filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, ".pdf", sep = "", 
                collapse = "")
            pdf(file = filename, height = height, width = width)
        }
    }
    
    cmap <- c("#AAAAFF", "#111166")
    
    # GSEA.HeatMapPlot2(V = A, row.names = A.row.names, col.names = A.names, main = 'Leading Subsets Assignment
    # (clustered)', sub = paste(doc.string, ' - ', phen2, sep=''), xlab=' ', ylab=' ', color.map = cmap)
    GSEA.HeatMapPlot2(V = t(A), row.names = A.names, col.names = A.row.names, main = "Leading Subsets Assignment (clustered)", 
        sub = paste(doc.string, " - ", phen2, sep = ""), xlab = " ", ylab = " ", color.map = cmap)
    
    text.filename <- paste(directory, doc.string, ".leading.assignment.clustered.", phen2, ".txt", sep = "", collapse = "")
    line.list <- c("Gene", A.row.names)
    line.header <- paste(line.list, collapse = "\t")
    line.length <- length(A.row.names) + 1
    write(line.header, file = text.filename, ncolumns = line.length)
    write.table(t(A), file = text.filename, append = T, quote = F, col.names = F, row.names = T, sep = "\t")
    
    if (non.interactive.run == F) {
        if (.Platform$OS.type == "windows") {
            savePlot(filename = filename, type = "jpeg", device = dev.cur())
            dev.off()
        } else if (.Platform$OS.type == "unix") {
            dev.off()
        }
    } else {
        dev.off()
    }
    
}
