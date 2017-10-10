#' plot overall/within group mean correlation coefficient of samples
#'
#' A function to calculate overall/within group mean correlation coefficient of samples for the Basic QC
#'
#'
#' @param expr  a expression matrix, samples in columns, expression values in rows:
#' @param meta a vector, for group, sample and color of group
#' @return correlation matrix and correlation figure
#' @export
#' @keywords correlation
#' @examples
#' #' cc_plot()


cc_plot <- function(expr, meta, ps = 8, fs = 23) {
    
    # sort by group and sample
    meta = meta[order(meta$Group, meta$Sample), ]
    expr = expr[, meta$Sample]
    grp = unique(meta$Group)
    smp = unique(meta$Sample)
    
    # within group Corr. Coefficient
    cc_grp = c()
    for (g in grp) {
        tmp = cor(as.matrix(expr[, meta[meta$Group == g, "Sample"]]), use = "pairwise")
        cc_grp = c(cc_grp, (apply(tmp, 2, sum, na.rm = T) - 1)/(apply(!is.na(tmp), 2, sum) - 1))
    }
    
    # overall Corr. Coefficient
    tmp = cor(as.matrix(expr[, smp]), use = "pairwise")
    cc_all = (apply(tmp, 2, sum, na.rm = T) - 1)/(apply(!is.na(tmp), 2, sum) - 1)
    
    # data frame for plotting
    dfp = meta
    dfp$sid = 1:nrow(dfp)
    dfp$cc_grp = cc_grp
    dfp$cc_all = cc_all
    dfp$Sample = factor(smp, levels = smp)
    
    plt = ggplot(dfp) + labs(title = "Correlation Plot: Within Group (red line) and Overall (blue line)") + geom_point(aes(x = Sample, 
        y = cc_grp, color = Group), size = ps) + geom_line(aes(x = sid, y = cc_grp), color = "red", size = 1) + 
        geom_point(aes(x = Sample, y = cc_all, color = Group), shape = 17, size = ps) + geom_line(aes(x = sid, y = cc_all), 
        color = "blue", size = 1) + scale_color_manual(values = unique(meta$Color)) + ylab("Corr. Coefficient") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = fs))
    # annotate('text', hjust=0, x=length(smp)/2, y=0.98, label='- Within Group', color='red') annotate('text',
    # hjust=0, x=length(smp)/2, y=0.95, label='- Overall', color='blue')
    
    return(list(plt = plt, dfp = dfp))
}

