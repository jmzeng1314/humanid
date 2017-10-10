#' write some files from the DEG results from topTable(limma)
#'
#'
#' @param DEG  DEG=topTable(fit,coef=2,adjust='BH')
#' @param imageType choose png or pdf or EMF
#' @param prefix A standard study ID
#' @return plot 3 kind of volcanic figures
#' @export
#' @keywords Volcanic_DEG
#' @examples
#' #'  Volcanic_DEG(DEG)
#'
#'
Volcanic_DEG <- function(DEG, imageType = "png", prefix = "test") {
    library(ggplot2)
    library(devEMF)
    res <- DEG
    res$Gene <- rownames(res)
    # Make a basic volcano plot
    logFC_Cutof <- with(res, mean(abs(logFC)) + 2 * sd(abs(logFC)))
    
    if (imageType == "png") {
        png(paste0(prefix, "_basic_volcanic.png"))
    } else if (imageType == "pdf") {
        pdf(paste0(prefix, "_basic_volcanic.pdf"))
    } else if (imageType == "EMF") {
        emf(paste0(prefix, "_basic_volcanic.EMF"))
    } else {
        stop("we just afford png ,emf or pdf for the volcanic figures")
    }
    
    with(res, plot(logFC, -log10(P.Value), pch = 20, main = "", xlab = "", ylab = ""))
    abline(h = -log10(0.05), col = "blue")
    abline(v = logFC_Cutof, col = "red")
    abline(v = -logFC_Cutof, col = "red")
    title(main = paste0("cutoff for logFC is ", round(logFC_Cutof, 3)), xlab = "log2 fold change", ylab = "-log10 p-value")
    dev.off()
    
    if (imageType == "png") {
        png(paste0(prefix, "_marker_gene_volcanic.png"))
    } else if (imageType == "pdf") {
        pdf(paste0(prefix, "_marker_gene_volcanic.pdf"))
    } else if (imageType == "EMF") {
        emf(paste0(prefix, "_marker_gene_volcanic.EMF"))
    } else {
        stop("we just afford png ,emf or pdf for the volcanic figures")
    }
    
    
    with(res, plot(logFC, -log10(P.Value), pch = 20, main = "Volcano plot"))
    # Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
    with(subset(res, adj.P.Val < 0.05), points(logFC, -log10(P.Value), pch = 20, col = "red"))
    with(subset(res, abs(logFC) > 1), points(logFC, -log10(P.Value), pch = 20, col = "orange"))
    with(subset(res, adj.P.Val < 0.05 & abs(logFC) > 1), points(logFC, -log10(P.Value), pch = 20, col = "green"))
    # Label points with the textxy function from the calibrate plot
    library(calibrate)
    with(subset(res, adj.P.Val < 0.05 & abs(logFC) > 1), textxy(logFC, -log10(P.Value), labs = Gene, cex = 0.8))
    dev.off()
    
    
    logFC_Cutof <- with(res, mean(abs(logFC)) + 2 * sd(abs(logFC)))
    res$change = as.factor(ifelse(res$P.Value < 0.05 & abs(res$logFC) > logFC_Cutof, ifelse(res$logFC > logFC_Cutof, 
        "UP", "DOWN"), "NOT"))
    this_tile <- paste0("Cutoff for logFC is ", round(logFC_Cutof, 3), "\nThe number of up gene is ", nrow(res[res$change == 
        "UP", ]), "\nThe number of down gene is ", nrow(res[res$change == "DOWN", ]))
    
    if (imageType == "png") {
        png(paste0(prefix, "_ggplot_volcanic.png"))
    } else if (imageType == "pdf") {
        pdf(paste0(prefix, "_ggplot_volcanic.pdf"))
    } else if (imageType == "EMF") {
        emf(paste0(prefix, "_ggplot_volcanic.EMF"))
    } else {
        stop("we just afford png ,emf or pdf for the volcanic figures")
    }
    
    
    ## Construct the plot object
    g = ggplot(data = res, aes(x = logFC, y = -log10(P.Value), color = change)) + geom_point(alpha = 0.4, size = 1.75) + 
        theme_set(theme_set(theme_bw(base_size = 20))) + xlab("log2 fold change") + ylab("-log10 p-value") + ggtitle(this_tile) + 
        theme(plot.title = element_text(size = 15, hjust = 0.5)) + scale_colour_manual(values = c("blue", "black", 
        "red"))  ## corresponding to the levels(res$change)
    print(g)
    dev.off()
}



