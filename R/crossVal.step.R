#' Cross validation for variable selection by stepwise models
#'
#' This function evaluates variable selection by stepwise models by repeated split sampling
#' @param pred data frame or matrix of predictor variables
#' @param resp data frame with response variables, rownames have to correspond
#'   to `pred`
#' @param var response variable (column name or index of `resp`)
#' @param rounds number of split sampling runs. Default `200`
#' @param frac fraction of samples used for training set. Default `0.7`
#' @param plot logical. plot observed vs predicted (as list of grobs to be
#'   plotted with [grid.arrange()]). Default `FALSE`
#' @param method method for correlation between observed and predicted, as in
#'   [cor()]. Default `pearson`
#' @param direction direction for stepwise selection as in [MASS::stepAIC]. `"forward"`, `"backward"` or `"both"`.
#' Default `"both"`
#' @param family family for glm as in [glm()]
#' @param ignoreError logical. If `TRUE` an iteration with error is skipped and a new subsample is taken; proceeded until desired number of repeats is reached; 
#' useful e.g for non-gaussian glms that do not converge for the training set. Default `FALSE`
#' @returns 
#' A list with five or six elements
#' * `$cor` correlations observed vs predicted 
#' * `$MSE` Mean squared errors observed vs predicted:
#' * `$summary` data frame with mean and sd of all cross validation runs
#' * `$plots` plots (list of grobs)
#' * `$variables` For each cross validation run: variables that appear in the model
#' * `$variablesFraction` fraction of cross validation runs each variable appears in the model
#' @details The function follows the same principle as [crossVal()], but includes stepwise variable selection in each iteration. This allows to assess the fraction
#' of cross validation runs a variable enters the model, and thus if the set of predictors depends on the dataset chosen. 
#' Plots of observed vs predicted values can be visualized with [gridExtra::grid.arrange()], e.g. `grid.arrange(grobs=object$plots)`
#' @seealso 
#' * [crossVal()] for cross validation of a final model
#' * [crossVal.ensemble()] for cross validation of a weighted ensemble model
#' @importFrom ggplot2 theme element_blank element_rect element_text ggplot aes geom_point geom_abline 
#' @importFrom egg set_panel_size
#' @importFrom grid unit
#' @importFrom randomForest randomForest
#' @importFrom stats glm predict cor sd
#' @importFrom MASS stepAIC
#' @export
#' @examples
#' set.seed(12)
#' cv <- crossVal.step(pred = MicEnvMod::pr, resp = MicEnvMod::resp,
#' var = "optimum_temperature", plot = TRUE, direction = "both")
#' # summary
#' cv$summary
#' # frequency variables appear in model
#' cv$variablesFraction
#' # plot
#' gridExtra::grid.arrange(grobs = cv$plots, ncol = 15)


crossVal.step <- function(pred, resp, var, rounds = 200, frac = 0.7, plot = FALSE, method = "pearson", direction = "both", family = "gaussian", ignoreError = FALSE){
  res <- list()
  res$cor <- c() # vector for correlation
  res$MSE <- c() # vector for means squared error
  Vars <- list()
  
  if(plot==TRUE){
    # plot settings
    plot.theme1 <- theme(panel.grid.major = element_blank(),
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
    res$plots <- list()
  }
  
  if (ignoreError == FALSE){
    for (i in 1:rounds){
      
      # repeat split sampling until all levels of categorical variables in test data are also in training
      repeat {
        ntrain <- round(0.7 * nrow(pred))
        sub <- sample(rownames(pred), replace = FALSE, size = ntrain)
        
        train_p <- pred[sub,]
        train_r <- resp[sub,var]
        test_p <- pred[!(rownames(pred)%in%sub),]
        test_r <- resp[!(rownames(pred)%in%sub), var]
        # check if for all categorical variables condition is met
        # names of categorical variables
        cat <- colnames(train_p)[-which(unlist(lapply(train_p, is.numeric), use.names = FALSE))]
        levTest <- c()
        for (j in cat){
          levTest[j] <- all(test_p[,j]%in%train_p[,j])
        }
        
        # exit if the condition is met for all categorical variables
        if (all(levTest)) break
      }
      
      # create model with training set
      # if direction is forward start with intercept only modell, if backward or both start with full model
      
      null <- glm(train_r ~ 1, data = train_p, family = family)
      full <- glm(train_r ~., data = train_p, family = family)
      if (direction == "forward"){
        mod <- stepAIC(null, scope = list(upper = full, lower = null), direction = direction, trace = FALSE)
      }else{
        mod <- stepAIC(full, scope = list(upper = full, lower = null), direction = direction, trace = FALSE)
      }
      
      # assess observed vs predicted in test set
      # variables selected in each run
      Vars[[i]] <- colnames(mod$model)[2:length(colnames(mod$model))]
      
      # data frame with observed and predicted values for test data set
      dat <- data.frame(Observed = test_r, Predicted = predict(mod, newdata = test_p, type = "response"))
      
      res$cor[i] <- cor(dat$Observed, dat$Predicted, method = method)
      res$MSE[i] <- mean((dat$Observed - dat$Predicted)^2)
      
      if(plot == TRUE){
        p <- ggplot(dat, aes(x = Observed, y = Predicted)) + geom_point() +
          geom_abline() + plot.theme1
        res$plots[[i]] <- set_panel_size(p, width = unit(1.5, "inch"), height = unit(1.5, "inch"))
      }
    }
  }else{
    i <- 1
    repeat{
      tryCatch(
        {
          # repeat split sampling until all levels of categorical variables in test data are also in training
          repeat {
            ntrain <- round(0.7 * nrow(pred))
            sub <- sample(rownames(pred), replace = F, size = ntrain)
            
            train_p <- pred[sub,]
            train_r <- resp[sub,var]
            test_p <- pred[!(rownames(pred)%in%sub),]
            test_r <- resp[!(rownames(pred)%in%sub), var]
            # check if for all categorical variables condition is met
            # names of categorical variables
            cat <- colnames(train_p)[-which(unlist(lapply(train_p, is.numeric), use.names = FALSE))]
            levTest <- c()
            for (j in cat){
              levTest[j] <- all(test_p[,j]%in%train_p[,j])
            }
            
            # exit if the condition is met for all categorical variables
            if (all(levTest)) break
          }
          
          # create model with training set
          # if direction is forward start with intercept only modell, if backward or both start with full model
          
          null <- glm(train_r ~ 1, data = train_p, family = family)
          full <- glm(train_r ~., data = train_p, family = family)
          if (direction == "forward"){
            mod <- stepAIC(null, scope = list(upper = full, lower = null), direction = direction, trace = FALSE)
          }else{
            mod <- stepAIC(full, scope = list(upper = full, lower = null), direction = direction, trace = FALSE)
          }
          
          # assess observed vs predicted in test set
          # variables selected in each run
          Vars[[i]] <- colnames(mod$model)[2:length(colnames(mod$model))]
          
          # data frame with observed and predicted values for test data set
          dat <- data.frame(Observed = test_r, Predicted = predict(mod, newdata = test_p, type = "response"))
          
          res$cor[i] <- cor(dat$Observed, dat$Predicted, method = method)
          res$MSE[i] <- mean((dat$Observed - dat$Predicted)^2)
          
          if(plot == TRUE){
            p <- ggplot(dat, aes(x = Observed, y = Predicted)) + geom_point() +
              geom_abline() + plot.theme1
            res$plots[[i]] <- set_panel_size(p, width = unit(1.5, "inch"), height = unit(1.5, "inch"))
          }
          
          i <- i+1
          
        }, error=function(e){
          cat("repeating with new data subset: (",conditionMessage(e), ")\n") 
        }
      )
      if (length(res$cor)==rounds) break
    }
  }
  
  suppressWarnings(Vars <- do.call(cbind, Vars))
  # different runs have different lengths and variables become recycled
  # warnings suppressed with suppressWarnings
  # replace recycled values with NA
  for(i in 1:ncol(Vars)){
    Vars[,i] <- c(Vars[!duplicated(Vars[,i]),i], rep(NA,length(which(duplicated(Vars[,i])))))
  }
  res$variables <- Vars
  # calculate fraction of runs where each predictor variable appears in model
  Vars <- as.vector(Vars)
  vec <- c()
  for (i in unique(Vars)[!is.na(unique(Vars))]){
    vec[i] <- length(which(Vars==i))/rounds
  }
  vec <- sort(vec, decreasing = T)
  
  res$variablesFraction <- vec
  
  res$summary <- c(meanCor = mean(res$cor), sdCor = sd(res$cor), meanMSE = mean(res$MSE), sdMSE = sd(res$MSE))
  return(res)
}