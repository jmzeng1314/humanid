#' Do differential expression analysis for the expression matrix
#'
#' please make sure that there are just two group for you expression matrix, you should modify the code in this function if there's more than 2 group.
#' I highly suggest you to use Limma to analysis,it's very easy, but you should make sure that
#' relevel(factor(group_list),'control')
#'
#' @param prefix    The prefix for all of the output files.( we don't need it now actually)
#' @param exprSet    Matrix with microarray expression values.
#' @param group_list Factors for two groups that are tested for differential expression.please make sure that relevel(factor(group_list),'control')
#' @param method     which method do you choose, limma,samr,rowttest in genefilter,resamplingPvalues in BioNet,t.test
#' @param destdir    where to store the files just download.( we don't need it now actually)
#' @return a data.frame just like topTable
#' @export
#' @keywords do_DEG_2groups
#' @examples
#' #' do_DEG_2groups('CLL',exprSet,group_list,method='limma')


do_DEG_2groups <- function(prefix = "GSE1009", exprSet = example_exprSet, group_list, method, destdir = ".") {
    # library(CLL) data(sCLLex) suppressMessages(library(limma)) exprSet = exprs(sCLLex) pdata=pData(sCLLex)
    # group_list = pdata$Disease do_DEG_2groups('CLL',exprSet,group_list)
    
    
    if (method == "limma") {
        library(limma)
        design = model.matrix(~factor(group_list))
        fit = lmFit(exprSet, design)
        fit = eBayes(fit)
        return(topTable(fit, coef = 2, n = Inf))
        
    } else {
        dat = exprSet
        group1 = which(group_list == levels(group_list)[1])
        group2 = which(group_list == levels(group_list)[2])
        dat1 = dat[, group1]
        dat2 = dat[, group2]
        dat = cbind(dat1, dat2)
        if (method == "t.test") {
            library(pi0)
            pvals = matrix.t.test(dat, 1, length(group1), length(group2))
        } else if (method == "samr") {
            library(samr)
            data = list(x = exprSet, y = as.numeric(as.factor(group_list)), geneid = as.character(1:nrow(exprSet)), 
                genenames = rownames(exprSet), logged2 = TRUE)
            samr.obj <- samr(data, resp.type = "Two class unpaired", nperms = 1000)
            pv = samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
            # delta.table <- samr.compute.delta.table(samr.obj)
            del <- 0.3
            # samr.plot(samr.obj, del=.3) siggenes.table<- samr.compute.siggenes.table(samr.obj, del, data, delta.table) ##
            # a list contains up gene list and down gene list
            pvals = pv
        } else if (method == "rowttest") {
            library(genefilter)
            library(impute)
            t.test <- rowttests(exprSet, fac = group_list)
            t.test[1:10, ]
            pvals = t.test$p.value
        } else if (method == "resamplingPvalues") {
            
            library(BioNet)
            results <- resamplingPvalues(exprSet, group_list, alternative = "two.sided")
            head(results$p.values, 10)
            pvals = results$p.values
        } else {
            stop("we just affort limma,samr,rowttest in genefilter,resamplingPvalues in BioNet,t.test to do DEG analysis!!!")
        }
        
        p.adj = p.adjust(pvals, method = "BH")
        avg_1 = rowMeans(dat1)
        avg_2 = rowMeans(dat2)
        FC = avg_2/avg_1
        results = cbind(avg_1, avg_2, FC, pvals, p.adj)
        colnames(results) = c("avg_1", "avg_2", "logFC", "P.Value", "adj.P.Val")  ## make sure return the same column name with limma
        rownames(results) = rownames(dat)
        results = as.data.frame(results)
        return(results[order(results$P.Value), ])
        
    }
}

# head(do_DEG_2groups('CLL',exprSet,group_list,method='limma'),10)
# head(do_DEG_2groups('CLL',exprSet,group_list,method='samr'),10)
# head(do_DEG_2groups('CLL',exprSet,group_list,method='rowttest'),10)
# head(do_DEG_2groups('CLL',exprSet,group_list,method='resamplingPvalues'),10)
# head(do_DEG_2groups('CLL',exprSet,group_list,method='t.test'),10)


