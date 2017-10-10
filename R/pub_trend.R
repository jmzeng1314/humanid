#' visualization for the trend of reseach area based on some keywords
#'
#' you have to install the packages:rentrez,ggplot2,reshape2 before using this function .
#' Then you can search any topic by the keywords you offer,such as '3C','4C','5C','Hi-C' or 'H3K4ME1','H3K4ME2' ,'H3K4ME3'.
#' It should be the same result when you search the keyword in pubmed : https://www.ncbi.nlm.nih.gov/pubmed/?term=H3K4ME1
#' Also you can narrow down the period for this search , such as from 2010 to 2016, as you will .
#' Thanks HaiTao (my friend) for the contribution .
#'
#' @param years A simple vector to narrow down the search area, default : 2010:2016
#' @param keywords A simple vector contains the keywords,such as '3C','4C','5C','Hi-C','ChIA???PET','3D chromosome'
#' @return a list contains 4 elements, two figures and the paper number in pubmed database.
#' @export
#' @keywords pub_trend
#' @examples
#' #' tmp <- pub_trend();tmp;tmp <- pub_trend(keywords=c('H3K4ME1','H3K4ME2' ,'H3K4ME3' ));
#'
#'
#'
pub_trend <- function(years = 2010:2016, keywords = c("3C", "4C", "5C", "Hi-C", "ChIA???PET", "3D chromosome")) {
    library(rentrez)
    library(ggplot2)
    library(reshape2)
    # papers_by_year function
    papers_by_year <- function(years, search_term) {
        return(sapply(years, function(y) entrez_search(db = "pubmed", term = search_term, mindate = y, maxdate = y, 
            retmax = 0)$count))
    }
    total_papers <- papers_by_year(years, "")
    total_papers
    ## https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
    trend_data <- sapply(keywords, function(t) papers_by_year(years, t))
    trend_props <- trend_data/total_papers
    
    trend_df <- melt(data.frame(years, trend_props), id.vars = "years")
    p1 <- ggplot(trend_df, aes(years, value, colour = variable))
    p1 <- p1 + geom_line(size = 1) + scale_y_log10("number of papers")
    
    trend_df <- melt(data.frame(years, trend_data), id.vars = "years")
    p2 <- ggplot(trend_df, aes(years, value, colour = variable))
    p2 <- p2 + geom_line(size = 1) + scale_y_log10("numberofpapers")
    
    return(list(trend_data = trend_data, trend_props = trend_props, p1 = p1, p2 = p2))
    
}
