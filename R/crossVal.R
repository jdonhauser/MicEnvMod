#' Cross validation by repeated split sampling
#'
#'This function evaluates model performance by repeated split sampling based on correlations
#'between observed and predicted values.
#' @param pred data frame or matrix of predictor variables
#' @param resp data frame with response variables, rownames have to correspond
#'   to `pred`
#' @param var response variable (column name or index of `resp`)
#' @param rounds number of split sampling runs. Default `200`
#' @param frac fraction of samples used for training set. Default `0.7`
#' @param plot logical. plot observed vs predicted (as list of grobs to be
#'   plotted with grid.arrange()). Default `FALSE`
#' @param method method for correlation between observed and predicted, as in
#'   [cor()]. Default `pearson`
#' @param type model type; "rf" for random forest or "glm"
#' @param family family for glm as in [glm()]. Default `"gaussian"`
#' @param ignoreError logical. If `TRUE` an iteration with error is skipped and a new subsample is taken; proceeded until desired number of repeats is reached; 
#' useful e.g for non-gaussian glms that do not converge for the training set. Default `FALSE`
#' @returns A list with three or four elements
#' * `$cor` correlations observed vs predicted 
#' * `$MSE` Mean squared errors observed vs predicted:
#' * `$summary` data frame with mean and sd of all cross validation runs
#' * `$plots` plots (list of grobs)
#' @details This function splits the data in a test and training set by randomly drawing samples without replacement. The function builds the model with the training set,
#' predicts the response variable for the test test. Measures of accuracy include correlation between observed and predicted values as well as mean squared errors. The procedure is repeated n times
#' and results are averaged. Measures of accuracy can be used as weights in an ensemble model.
#' Plots of observed vs predicted values can be visualized with [gridExtra::grid.arrange()], e.g. `grid.arrange(grobs=object$plots)`
#' @references Thuiller, W., Lafourcade, B., Engler, R., & Araujo, M. B. (2009). 
#' BIOMOD - A platform for ensemble forecasting of species distributions.
#' *Ecography*, 32(3), 369-373. https://doi.org/10.1111/j.1600-0587.2008.05742.x
#' @seealso 
#' * [crossVal.step()] for cross validation with stepwise selection
#' * [crossVal.ensemble()] for cross validation of a weighted ensemble model
#' @importFrom ggplot2 theme element_blank element_rect element_text ggplot aes geom_point geom_abline 
#' @importFrom egg set_panel_size
#' @importFrom grid unit
#' @importFrom randomForest randomForest
#' @importFrom stats glm predict cor sd
#' @export
#' @examples
#' # for random forest
#' set.seed(12)
#' rf <- randomForest::randomForest(x = MicEnvMod::pr, y = resp[,2],
#' ntree = 1000, importance = TRUE, nPerm = 100)
#' cv <- crossVal(pred = MicEnvMod::pr, resp = MicEnvMod::resp,
#' var = "optimum_temperature", plot = TRUE, type = "rf")
#' # summary
#' cv$summary
#' # plot
#' gridExtra::grid.arrange(grobs = cv$plots, ncol = 15)
#' 
#' # for stepwise glm
#' full <- glm(MicEnvMod::resp[,2] ~., data = MicEnvMod::pr, family = gaussian)
#' null <- glm(MicEnvMod::resp[,2] ~1, data = MicEnvMod::pr, family = gaussian)
#' glmStep <- MASS::stepAIC(full, scope = list(upper = full, lower = null), direction = 'both')
#' set.seed(12)
#' cv <- crossVal(pred = glmStep$model[,2:ncol(glmStep$model)], resp = MicEnvMod::resp,
#' var = "optimum_temperature", plot = TRUE, type = "glm")
#' # summary
#' cv$summary
#' # plot
#' gridExtra::grid.arrange(grobs = cv$plots, ncol = 15)
 
crossVal <- function(pred, resp, var, rounds = 200, frac = 0.7, plot = FALSE, method = "pearson", type = "rf", family = "gaussian", ignoreError = FALSE){
  res <- list()
  res$cor <- c() # vector for correlation
  res$MSE <- c() # vector for means squared error
  
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
  
  if (ignoreError==FALSE){
    for (i in 1:rounds){
      repeat {
        ntrain <- round(0.7 * nrow(pred))
        sub <- sample(rownames(pred), replace = FALSE, size = ntrain)
        
        train_p <- pred[sub,, drop = FALSE]
        train_r <- resp[sub,var]
        test_p <- pred[!(rownames(pred)%in%sub),,drop = FALSE]
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
      
      
      if (type == "rf"){
        # train model on training set
        mod <- randomForest(x = train_p, y = train_r, ntree = 1000, importance = FALSE)
        # data frame with observed and predicted values for test data set
        dat <- data.frame(Observed = test_r, Predicted = predict(mod, newdata = test_p))
      }
      
      if (type == "glm"){
        
        # create model with training set
        mod <- glm(train_r ~., data = train_p, family = family)
        # data frame with observed and predicted values for test data set
        dat <- data.frame(Observed = test_r, Predicted = predict(mod, newdata = test_p, type = "response"))
      }
      
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
          repeat {
            ntrain <- round(0.7 * nrow(pred))
            sub <- sample(rownames(pred), replace = FALSE, size = ntrain)
            
            train_p <- pred[sub,, drop = FALSE]
            train_r <- resp[sub,var]
            test_p <- pred[!(rownames(pred)%in%sub),,drop = FALSE]
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
          
          
          if (type == "rf"){
            # train model on training set
            mod <- randomForest(x = train_p, y = train_r, ntree = 1000, importance = FALSE)
            # data frame with observed and predicted values for test data set
            dat <- data.frame(Observed = test_r, Predicted = predict(mod, newdata = test_p))
          }
          
          if (type == "glm"){
            
            # create model with training set
            mod <- glm(train_r ~., data = train_p, family = family)
            # data frame with observed and predicted values for test data set
            dat <- data.frame(Observed = test_r, Predicted = predict(mod, newdata = test_p, type = "response"))
          }
          
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
  
  res$summary <- c(meanCor = mean(res$cor), sdCor = sd(res$cor), meanMSE = mean(res$MSE), sdMSE = sd(res$MSE))
  return(res)
}
