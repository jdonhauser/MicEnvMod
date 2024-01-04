
#' Plot variable importance for random forest with ggplot
#'
#' @param mod a model of class `randomForest`
#' @param pt optional settings for plots. a `ggplot` `theme`
#' @param decreasing logical. sort variable importance in decreasing order. Default = `FALSE`
#' @returns grob to be plotted with [gridExtra::grid.arrange()]
#' @importFrom ggplot2 theme element_blank element_rect element_text ggplot aes geom_col xlab
#' @importFrom egg set_panel_size
#' @importFrom grid unit
#' @importFrom randomForest varImpPlot
#' @export
#' @examples
#' # random forest model
#' set.seed(12)
#' rf <- randomForest::randomForest(x = MicEnvMod::pr, y = resp[,2], ntree = 1000,
#' importance = TRUE, nPerm = 100)
#' 
#' # Variance importance plot
#' p <- varImpPlotgg(rf, decreasing = TRUE)
#' gridExtra::grid.arrange(p)

varImpPlotgg <- function(mod, pt = NULL, decreasing = FALSE){
  if(is.null(pt)){
    # default plot theme
    pt <- theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_rect(fill = "white",
                                                         colour = "black",
                                                         size = 0.5, linetype = "solid"),
                         panel.border= element_rect(fill=NA,size = 0.5, linetype = 'solid',colour = "black"),
                         axis.text.x = element_text(size=13),axis.text.y = element_text(size=13),legend.text = element_text(size=13),
                         axis.title = element_text(size=14),
                         legend.title = element_text(color = "black", size = 14),
                         strip.text.x = element_text(size=14),
                         strip.background = element_rect(colour="black", fill="white")
    )
  }
  
  
  varImp <- varImpPlot(mod,sort = TRUE)
  varImp <- as.data.frame(varImp)
  varImp <- varImp[order(varImp$`%IncMSE`, decreasing = decreasing),]
  varImp$Variable <- factor(rownames(varImp),levels = rownames(varImp))
  p <- ggplot(varImp, aes(x = `%IncMSE`, y = Variable)) + geom_col() + 
    xlab('% increase MSE') + pt
  p <- set_panel_size(p, width = unit(2.5, "inch"), height = unit(nrow(varImp)/2, "inch"))
}
