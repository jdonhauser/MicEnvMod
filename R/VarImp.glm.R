# does not work if there is only one variable in model

#' Variable importance for a glm
#'
#' @param mod model of class `glm`
#' @param perm number of permutations
#' @param plot logical. create plot for variable importance. (to be
#'   plotted with [gridExtra::grid.arrange()]). Default `TRUE`
#' @param pt optional. settings for plots. a ggplot `theme`
#' @param decreasing logical. sort variable importance in decreasing order. Default = `FALSE`
#' @returns A list with one or two elements
#' * `$varImpDat` a data frame with one row for each predictor and two columns (mean and sd variable importance)
#' * `$res$varImpPlot` a grob, returned if plot = TRUE 
#' @details This function calculates variable importance by permuting the predictor variables one at a time and calculating
#' the % increase in mean squared error. The function is analogous to the variable importance created by [randomForest()].
#' Plots can be visualized with [grid.arrange()], e.g. `grid.arrange(object$varImpPlot)`
#' @importFrom ggplot2 theme element_blank element_rect element_text ggplot aes geom_col xlab
#' @importFrom egg set_panel_size
#' @importFrom grid unit
#' @importFrom stats glm sd
#' @importFrom permute shuffle
#' @export
#' @examples
#' set.seed(12)
#' # stepwise glm
#' full <- glm(MicEnvMod::resp[,2] ~., data = MicEnvMod::pr, family = gaussian)
#' null <- glm(MicEnvMod::resp[,2] ~1, data = MicEnvMod::pr, family = gaussian)
#' glmStep <- MASS::stepAIC(full, scope = list(upper = full, lower = null),
#' direction = 'both')
#' 
#' # Variance importance 
#' varImp <- VarImp.glm(mod = glmStep, perm = 100, plot = TRUE)
#' # plot
#' gridExtra::grid.arrange(varImp$varImpPlot)

VarImp.glm <- function(mod, perm = 100,plot = TRUE, pt = NULL, decreasing = FALSE){
  
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
  
  # predictors in final model
  pr <- colnames(mod$model)[2:length(colnames(mod$model))]
  # mean squared error original model
  mse <- mean((mod$fitted.values - mod$y)^2)
  
  varImp <- c()
  varImpSD <- c()
  
  for (i in pr){
    tmpVar <- c()
    for(j in 1:perm){
      dat <- mod$model[,2:length(colnames(mod$model)),drop = FALSE]
      dat[,i] <- dat[,i][shuffle(dat[,i])]
      m <- glm(mod$y ~., data = dat, family = mod$family$family)
      mse_i <- mean((m$fitted.values - mod$y)^2)
      # % increase in mean squared error when variable i is permuted
      tmpVar[j] <- ((mse_i-mse)/mse)*100
    }
    # mean of repetitions
    varImp[i] <- mean(tmpVar)
    varImpSD[i] <- sd(tmpVar)
  }
  
  varImp <- data.frame(Importance = varImp, ImportanceSD = varImpSD)
  varImp <- varImp[order(varImp$Importance, decreasing = decreasing),]
  varImp$Variable <- factor(rownames(varImp),levels = rownames(varImp))
  
  res <- list()
  res$varImpDat <- varImp
  
  if (plot == TRUE){
    # res$varImpPlot <- list()
    p <- ggplot(varImp, aes(x = Importance, y = Variable)) + geom_col() + 
      xlab('% increase MSE') + pt
    p <- set_panel_size(p, width = unit(2.5, "inch"), height = unit(nrow(varImp)/2, "inch"))
    
    res$varImpPlot <- p
  }
  return(res)
}