#' Monovariate response plots
#'
#' This function creates monovariate response plots to evaluate how a variable
#' responds to a particular predictor
#' @param mod model for which to generate the response plot.
#' @param pred table of predictor variables
#' @param len number of data points for which response is predicted. Default `200`
#' @param pt optional graphic settings. A ggplot theme can be provided
#' @param w width of the plot panel (inch) 
#' @param h width of the plot panel (inch)
#' @param type full model or stepwise selected model (default `full`)
#' @param scale `"free"` or `"fixed"`. if set to `"free"` (default), the range of the y axis
#'   (response variable) varies across plots. if set to `"fixed"` the same scale is
#'   used for all plots
#' @param cols optional named color vector for bars representing levels of a
#'   categorical variable
#' @returns list of grobs, scatter plots for continuous variables and barplots for categorical variables as function of each predictor variable
#' @details
#' This function creates response plots by setting all predictor variables to their mean (for continuous variable) or the first level, except the variable of interest and 
#' computes predictions for the response variable. 
#' Plots can be visualized with [gridExtra::grid.arrange()], e.g. `grid.arrange(grobs=object)`
#' @references Elith, J., Ferrier, S., Huettmann, F., & Leathwick, J. (2005). The evaluation strip: 
#' A new and robust method for plotting predicted responses from species distribution models. *Ecological Modelling*,
#'  186(3), 280-289. https://doi.org/10.1016/j.ecolmodel.2004.12.007
#' @importFrom ggplot2 theme element_blank element_rect element_text ggplot aes geom_line scale_y_continuous xlab
#' geom_col scale_fill_manual 
#' @importFrom egg set_panel_size
#' @importFrom grid unit
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' # for random forest
#' set.seed(12)
#' rf <- randomForest::randomForest(x = MicEnvMod::pr, y = resp[,2], ntree = 1000,
#'  importance = TRUE, nPerm = 100)
#' respM <- respMono(mod = rf, pred = MicEnvMod::pr, scale = "fixed")
#' gridExtra::grid.arrange(grobs = respM, ncol = 3)
#' 
#' # for stepwise glm
#' full <- glm(MicEnvMod::resp[,2] ~., data = MicEnvMod::pr, family = gaussian)
#' null <- glm(MicEnvMod::resp[,2] ~1, data = MicEnvMod::pr, family = gaussian)
#' glmStep <- MASS::stepAIC(full, scope = list(upper = full, lower = null),
#' direction = 'both')
#' respM <- respMono(mod = glmStep, pred = MicEnvMod::pr, scale = "fixed")
#' gridExtra::grid.arrange(grobs = respM, ncol = 3)

respMono <- function(mod, pred, len = 200, pt = NULL, w = 2, h = 2, type = "full", scale = "free", cols = NULL){
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
  resp <- list() # list for storing plot data
  gl <- list() # list for storing plots
  
  # vector of numeric columns (for factor plot would have to be different, could be implemented later)
  # vars <- colnames(pred[,which(unlist(lapply(pred, is.numeric), use.names = FALSE))])
  
  vars <- colnames(pred)
  
  if (type == "step"){
    pr <- colnames(mod$model)[2:length(colnames(mod$model))]
    vars <- vars[vars%in%pr]
  }
  
  if (scale == "free"){
    for (k in vars){
      # set all continuous variables to mean and factors to first value 
      # create vector with means or first value for all variables except variable of interest
      m <- c()
      for (i in colnames(pred)[-which(colnames(pred)==k)]){
        if(is.numeric(pred[,i])){
          m[i] <- mean(pred[,i])
        }else{
          m[i] <- as.character(unique(pred[,i])[1])
        }
      }
      
      # create table with len values where for variable of interest range is used and for rest mean
      # start with variable of interst
      # for numeric create len points along range, for factor number of levels
      if(is.numeric(pred[,k])){
        nd <- data.frame(k = seq(min(pred[,k]),max(pred[,k]),length.out = len))
        colnames(nd)[1] <- k
      }else{
        nd <- data.frame(k = unique(pred[,k]))
        colnames(nd)[1] <- k
      }
      
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
      
      if(is.numeric(pred[,k])){
        p <- ggplot(nd, aes(x = get(k), y = prediction)) + geom_line() +
          xlab(k) +pt
        p <- set_panel_size(p, width = unit(w, "inch"), height = unit(h, "inch"))
        gl[[k]] <- p
      }else{
        # sort factor by prediction
        nd[,k] <- factor(nd[,k], levels = nd[,k][order(nd[,'prediction'])])
        # color for factor levels
        if(is.null(cols)){
          values = rep("grey30",length(levels(nd[,k])))
        }else{
          values = cols[levels(nd[,k])]
        }
        
        p <- ggplot(nd, aes(x = get(k), y = prediction, fill = get(k))) + geom_col(show.legend = FALSE) +
          scale_fill_manual(values = values) +
          xlab(k) + pt + theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
        p <- set_panel_size(p, width = unit(w, "inch"), height = unit(h, "inch"))
        gl[[k]] <- p
      }
    }
  }else{
    
    # determine minimum and maximum values for prediction for all plots
    minvals <- c()
    maxvals <- c()
    
    for (k in vars){
      # set all continuous variables to mean and factors to first value 
      # create vector with means or first value for all variables except both  variables of interest
      
      m <- c()
      for (i in colnames(pred)[-which(colnames(pred)==k)]){
        if(is.numeric(pred[,i])){
          m[i] <- mean(pred[,i])
        }else{
          m[i] <- as.character(unique(pred[,i])[1])
        }
      }
      
      # create table with len values where for variable of interest range is used and for rest mean
      # start with variable of interst
      # for numeric create len points along range, for factor number of levels
      if(is.numeric(pred[,k])){
        nd <- data.frame(k = seq(min(pred[,k]),max(pred[,k]),length.out = len))
        colnames(nd)[1] <- k
      }else{
        nd <- data.frame(k = unique(pred[,k]))
        colnames(nd)[1] <- k
      }
      
      
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
      resp[[k]] <- nd
      
      minvals[k] <- min(nd$prediction)
      maxvals[k] <- max(nd$prediction)
      
    } 
    
    # plot with min and max for prediction derived in previous loop
    for (k in vars){
      
      if(is.numeric(pred[,k])){
        p <- ggplot(resp[[k]], aes(x = get(k), y = prediction)) + geom_line() +
          scale_y_continuous(limits = c(min(minvals),max(maxvals))) + 
          xlab(k) + pt 
        p <- set_panel_size(p, width = unit(w, "inch"), height = unit(h, "inch"))
        gl[[k]] <- p
      }else{
        # sort factor by prediction
        resp[[k]][,k] <- factor(resp[[k]][,k], levels = resp[[k]][,k][order(resp[[k]][,'prediction'])])
        # color for factor levels
        if(is.null(cols)){
          values = rep("grey30",length(levels(resp[[k]][,k])))
        }else{
          values = cols[levels(resp[[k]][,k])]
        }
        
        p <- ggplot(resp[[k]], aes(x = get(k), y = prediction, fill = get(k))) + geom_col(show.legend = FALSE) +
          scale_fill_manual(values = values) +
          coord_cartesian(ylim= c(min(minvals),max(maxvals))) +
          xlab(k) + pt + theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
        p <- set_panel_size(p, width = unit(w, "inch"), height = unit(h, "inch"))
        gl[[k]] <- p
        
      }
    } 
  }
  
  
  return(gl)
}