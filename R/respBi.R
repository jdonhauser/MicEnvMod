#' Bivariate response plots
#'
#' This function creates bivariate response plots to evaluate how a variable
#' responds to two predictors at the same time
#' @param mod model for which to generate the response plot.
#' @param pred table of predictor variables
#' @param len number of data points for which response is predicted. Default `50`
#' @param pt graphic settings. A ggplot theme can be provided
#' @param w width of the plot panel (inch) 
#' @param h width of the plot panel (inch)
#' @param type full model or stepwise selected model (default `full`)
#' @returns list of grobs. heatmaps for response variable as a function of all pairwise combinations of predictors
#' @details
#' This function creates response plots by setting all predictor variables to their mean except the two variables of interest and 
#' computes predictions for the response variable.
#' Plots can be visualized with [gridExtra::grid.arrange()], e.g. `grid.arrange(grobs=object)`
#' @references Elith, J., Ferrier, S., Huettmann, F., & Leathwick, J. (2005). The evaluation strip: 
#' A new and robust method for plotting predicted responses from species distribution models. *Ecological Modelling*,
#'  186(3), 280-289. https://doi.org/10.1016/j.ecolmodel.2004.12.007
#' @importFrom ggplot2 theme element_blank element_rect element_text ggplot aes geom_raster xlab ylab
#' @importFrom egg set_panel_size
#' @importFrom grid unit
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @importFrom viridis scale_fill_viridis
#' @export
#' @seealso [respMono()]
#' @examples
#' # for random forest
#' set.seed(12)
#' rf <- randomForest::randomForest(x = MicEnvMod::pr, y = resp[,2], ntree = 1000,
#'  importance = TRUE, nPerm = 100)
#' respB <- respBi(mod = rf, pred = MicEnvMod::pr)
#' gridExtra::grid.arrange(grobs = respB, ncol = 3)
#'
#' # for stepwise glm
#' full <- glm(MicEnvMod::resp[,2] ~., data = MicEnvMod::pr, family = gaussian)
#' null <- glm(MicEnvMod::resp[,2] ~1, data = MicEnvMod::pr, family = gaussian)
#' glmStep <- MASS::stepAIC(full, scope = list(upper = full, lower = null),
#' direction = 'both')
#' respM <- respMono(mod = glmStep, pred = MicEnvMod::pr, scale = "fixed")
#' gridExtra::grid.arrange(grobs = respM, ncol = 7)


respBi <- function(mod, pred, len = 50, pt = NULL, w = 3, h = 3, type = "full"){
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
  
  gl <- list() # list for storing plots 
  
  # vector of numeric columns (for factor plot would have to be different, could be implemented later)
  vars <- colnames(pred[,which(unlist(lapply(pred, is.numeric), use.names = FALSE))])
  
  if (type == "step"){
    pr <- colnames(mod$model)[2:length(colnames(mod$model))]
    vars <- vars[vars%in%pr]
  }
  
  for (k in vars){
    for(l in vars){
      # set all continuous variables to mean and factors to first value 
      # create vector with means or first value for all variables except both  variables of interest
      m <- c()
      for (i in colnames(pred)[c(-which(colnames(pred)==k),-which(colnames(pred)==l))]){
        # for (i in colnames(pred)[c(-grep(k,colnames(pred)),-grep(l,colnames(pred)))]){
        if(is.numeric(pred[,i])){
          m[i] <- mean(pred[,i])
        }else{
          m[i] <- as.character(unique(pred[,i])[1])
        }
      }
      
      # create all combinations of values of the two variables
      nd <- data.frame(k = rep(seq(min(pred[,k]),max(pred[,k]),length.out = len), len),
                       l = rep(seq(min(pred[,l]),max(pred[,l]),length.out = len), each = len))
      
      colnames(nd) <- c(k,l)
      
      # add other variables
      # convert numberic variables back to numeric, rest becomes / remains character with rep
      for (i in names(m)){
        if(is.numeric(pred[,i])){
          nd[,i] <- as.numeric(rep(m[i],nrow(nd)))
        }else{
          nd[,i] <- rep(m[i],nrow(nd))
          nd[,i] <- factor(nd[,i], levels = levels(pred[,i]))
        }
      }
      
      nd$prediction <- predict(mod, newdata=nd, type = "response")
      
      p <- ggplot(nd, aes(x = get(k), y = get(l), fill = prediction)) + geom_raster() +
        xlab(k) + ylab(l)+
        scale_fill_viridis(limits = c(min(nd$prediction), max(nd$prediction))) +
        pt
      p <- set_panel_size(p, width = unit(w, "inch"), height = unit(h, "inch"))
      gl[[paste(k,l,sep = '_')]] <- p
      
    }
  }

  
  return(gl)
}
