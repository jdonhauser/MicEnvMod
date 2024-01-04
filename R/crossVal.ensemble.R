#' Cross validation by repeated split sampling for ensemble models
#'
#'This function evaluates performance of a weighted ensemble model by repeated split sampling 
#'based on correlations between observed and predicted values.
#' @param mods named list of all models to be included in the ensemble. Currently models of class `glm` and 
#' `randomForest` supported
#' @param resp data frame with response variables, rownames have to correspond
#'   to `pred`
#' @param var response variable (column name or index of `resp`)
#' @param rounds number of split sampling runs. Default `200`
#' @param frac fraction of samples used for training set. Default `0.7`
#' @param weights vector of weights for each model, has to be in same order as `mods`.
#' @param plot logical. plot observed vs predicted (as list of grobs to be
#'   plotted with [gridExtra::grid.arrange()]). Default `FALSE`
#' @param method method for correlation between observed and predicted, as in
#'   [cor()]. Default `pearson`
#' @param ignoreError logical. If `TRUE` an iteration with error is skipped and a new subsample is taken; proceeded until desired number of repeats is reached; 
#' useful e.g for non-gaussian glms that do not converge for the training set. Default `FALSE`
#' @returns A list with three or four elements
#' * `$cor` correlations observed vs predicted 
#' * `$MSE` Mean squared errors observed vs predicted:
#' * `$summary` data frame with mean and sd of all cross validation runs
#' * `$plots` plots (list of grobs)
#' @details Cross validation for a weighted ensemble model. This function creates splits the data in a test and training set by randomly drawing samples without replacement. The function all models with the training set,
#' predicts the response variable for the test set as a weighted average of predictions from the individual models. Accuracy measures from cross validation of the individual models can
#' be used as weight. Measures of accuracy include correlation between observed and predicted values as well as mean squared errors. The procedure is repeated n times
#' and results are averaged. 
#' Plots of observed vs predicted values can be visualized with [gridExtra::grid.arrange()], e.g. `grid.arrange(grobs=object$plots)`
#' @references Araujo, M. B., & New, M. (2007). Ensemble forecasting of species distributions. *Trends in Ecology and Evolution*, 
#' 22(1), 42-47. https://doi.org/10.1016/j.tree.2006.09.010
#' @seealso 
#' * [crossVal()] for cross validation of a single mode
#' * [crossVal.step()] for cross validation with stepwise selection
#' @importFrom ggplot2 theme element_blank element_rect element_text ggplot aes geom_point geom_abline 
#' @importFrom egg set_panel_size
#' @importFrom grid unit
#' @importFrom randomForest randomForest
#' @importFrom stats glm predict cor weighted.mean sd
#' @export
#' @examples
#' # cross validation for individual models to determine weights
#' pred <- MicEnvMod::pr
#' # model1: random forest
#' set.seed(12)
#' rf <- randomForest::randomForest(x = pred, y = resp[,2], ntree = 1000,
#' importance = TRUE, nPerm = 100)
#' cv <- crossVal(pred = pred, resp = MicEnvMod::resp,
#' var = "optimum_temperature", plot = FALSE, type = "rf")
#' # use correlation as weight
#' weights <- c()
#' weights["rf"] <- cv$summary[1]
#' 
#' # model2:  stepwise glm
#' full <- glm(MicEnvMod::resp[,2] ~., data = pred, family = gaussian)
#' null <- glm(MicEnvMod::resp[,2] ~1, data = pred, family = gaussian)
#' glmStep <- MASS::stepAIC(full, scope = list(upper = full, lower = null),
#' direction = 'both')
#' cv <- crossVal(pred = glmStep$model[,2:ncol(glmStep$model)], resp = MicEnvMod::resp,
#' var = "optimum_temperature", plot = FALSE, type = "glm")
#' weights["glm"] <- cv$summary[1]
#' 
#' # cross validation ensemble
#' # list of models
#' mods <- list(rf,glmStep)
#' names(mods) <- c("rf","glm")
#' 
#' cv <- crossVal.ensemble(mods = mods, resp = MicEnvMod::resp, var = "optimum_temperature",
#' plot = TRUE, weights = weights)
#' # summary
#' cv$summary
#' # plot
#' gridExtra::grid.arrange(grobs = cv$plots, ncol = 15)

crossVal.ensemble <- function(mods, resp, var, rounds = 200, frac = 0.7, weights, plot = FALSE, method = "pearson", ignoreError = F){
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
  
  if (ignoreError == FALSE){
    for (i in 1:rounds){
      train_p <- list()
      train_r <- list()
      test_p <- list()
      test_r <- list()
      
      # repeat split sampling until all levels of categorical variables in test data are also in training
      # split sampling subset defined for first model; for further models same subset used
      
      repeat {
        for (j in 1:length(mods)){
          if(j==1){
            if (all(class(mods[[j]]) == "randomForest" )){
              ntrain <- round(0.7 * nrow(get(mods[[j]]$call$x)))
              sub <- sample(rownames(get(mods[[j]]$call$x)), replace = FALSE, size = ntrain)
              
              train_p[[names(mods)[j]]] <- get(mods[[j]]$call$x)[sub,]
              train_r[[names(mods)[j]]] <- resp[sub,var]
              test_p[[names(mods)[j]]] <- get(mods[[j]]$call$x)[!(rownames(get(mods[[j]]$call$x))%in%sub),]
              test_r[[names(mods)[j]]] <- resp[!(rownames(get(mods[[j]]$call$x))%in%sub), var]
            }else{
              if ("glm"%in%class(mods[[j]])){
                ntrain <- round(0.7 * nrow(mods[[j]]$model[,2:ncol(mods[[j]]$model)]))
                sub <- sample(rownames(mods[[j]]$model[,2:ncol(mods[[j]]$model)]), replace = F, size = ntrain)
                
                train_p[[names(mods)[j]]] <- mods[[j]]$model[sub,2:ncol(mods[[j]]$model),drop=FALSE]
                train_r[[names(mods)[j]]] <- resp[sub,var]
                test_p[[names(mods)[j]]] <- mods[[j]]$model[!(rownames(mods[[j]]$model)%in%sub),2:ncol(mods[[j]]$model),drop=FALSE]
                test_r[[names(mods)[j]]] <- resp[!(rownames(mods[[j]]$model)%in%sub), var]
              }
            }
          }else{
            if (all(class(mods[[j]]) == "randomForest" )){
              train_p[[names(mods)[j]]] <- get(mods[[j]]$call$x)[sub,]
              train_r[[names(mods)[j]]] <- resp[sub,var]
              test_p[[names(mods)[j]]] <- get(mods[[j]]$call$x)[!(rownames(get(mods[[j]]$call$x))%in%sub),]
              test_r[[names(mods)[j]]] <- resp[!(rownames(get(mods[[j]]$call$x))%in%sub), var]
            }else{
              if ("glm"%in%class(mods[[j]])){
                train_p[[names(mods)[j]]] <- mods[[j]]$model[sub,2:ncol(mods[[j]]$model),drop=FALSE]
                train_r[[names(mods)[j]]] <- resp[sub,var]
                test_p[[names(mods)[j]]] <- mods[[j]]$model[!(rownames(mods[[j]]$model)%in%sub),2:ncol(mods[[j]]$model),drop=FALSE]
                test_r[[names(mods)[j]]] <- resp[!(rownames(mods[[j]]$model)%in%sub), var]
              }
            }
          }
        }
        
        
        # check if for all categorical variables condition is met (compare training vs test dat from all models)
        train_comb <- do.call(cbind, train_p)
        test_comb <- do.call(cbind, test_p)
        
        # names of categorical variables
        catVar <- colnames(train_comb)[-which(unlist(lapply(train_comb, is.numeric), use.names = FALSE))]
        levTest <- c()
        for (k in catVar){
          levTest[k] <- all(test_comb[,k]%in%train_comb[,k])
        }
        
        # exit if the condition is met for all categorical variables
        if (all(levTest)) break
      }
      
      # train models on training set
      modsNew <- list()
      for (l in names(mods)){
        if (all(class(mods[[l]]) == "randomForest" )){
          modsNew[[l]] <- randomForest(x = train_p[[l]], y = train_r[[l]], ntree = 1000, importance = FALSE)
        }else{
          if("glm"%in%class(mods[[l]])){
            modsNew[[l]] <- glm(train_r[[l]] ~., data = train_p[[l]], family = mods[[l]]$family$family)
          }
        }
      }
      # data frame with observed  for test data set and predicted values from each model
      dat <- data.frame(Observed = test_r[[1]])
      for(m in names(modsNew)){
        if("glm"%in%class(modsNew[[m]])){
          dat[,m] <- predict(modsNew[[m]], newdata = test_p[[m]], type = "response")
        }else{
          dat[,m] <- predict(modsNew[[m]], newdata = test_p[[m]])
        }
      }
      
      # calculate weighted ensemble prediction from individual prediction
      dat[,'ensemblePrediction'] <- rep(10,nrow(dat))
      for (n in 1:nrow(dat)){
        dat[n,'ensemblePrediction'] <- weighted.mean(dat[n,2:(ncol(dat)-1)],weights)
      }
      
      res$cor[i] <- cor(dat$Observed, dat$ensemblePrediction, method = method)
      res$MSE[i] <- mean((dat$Observed - dat$ensemblePrediction)^2)
      
      if(plot == T){
        p <- ggplot(dat, aes(x = Observed, y = ensemblePrediction)) + geom_point() +
          geom_abline() + plot.theme1
        res$plots[[i]] <- set_panel_size(p, width = unit(1.5, "inch"), height = unit(1.5, "inch"))
      }
    }
  }else{
    i <- 1
    repeat{
      tryCatch(
        {
          train_p <- list()
          train_r <- list()
          test_p <- list()
          test_r <- list()
          
          # repeat split sampling until all levels of categorical variables in test data are also in training
          # split sampling subset defined for first model; for further models same subset used
          
          repeat {
            for (j in 1:length(mods)){
              if(j==1){
                if (all(class(mods[[j]]) == "randomForest" )){
                  ntrain <- round(0.7 * nrow(get(mods[[j]]$call$x)))
                  sub <- sample(rownames(get(mods[[j]]$call$x)), replace = FALSE, size = ntrain)
                  
                  train_p[[names(mods)[j]]] <- get(mods[[j]]$call$x)[sub,]
                  train_r[[names(mods)[j]]] <- resp[sub,var]
                  test_p[[names(mods)[j]]] <- get(mods[[j]]$call$x)[!(rownames(get(mods[[j]]$call$x))%in%sub),]
                  test_r[[names(mods)[j]]] <- resp[!(rownames(get(mods[[j]]$call$x))%in%sub), var]
                }else{
                  if ("glm"%in%class(mods[[j]])){
                    ntrain <- round(0.7 * nrow(mods[[j]]$model[,2:ncol(mods[[j]]$model)]))
                    sub <- sample(rownames(mods[[j]]$model[,2:ncol(mods[[j]]$model)]), replace = FALSE, size = ntrain)
                    
                    train_p[[names(mods)[j]]] <- mods[[j]]$model[sub,2:ncol(mods[[j]]$model),drop=FALSE]
                    train_r[[names(mods)[j]]] <- resp[sub,var]
                    test_p[[names(mods)[j]]] <- mods[[j]]$model[!(rownames(mods[[j]]$model)%in%sub),2:ncol(mods[[j]]$model),drop=FALSE]
                    test_r[[names(mods)[j]]] <- resp[!(rownames(mods[[j]]$model)%in%sub), var]
                  }
                }
              }else{
                if (all(class(mods[[j]]) == "randomForest" )){
                  train_p[[names(mods)[j]]] <- get(mods[[j]]$call$x)[sub,]
                  train_r[[names(mods)[j]]] <- resp[sub,var]
                  test_p[[names(mods)[j]]] <- get(mods[[j]]$call$x)[!(rownames(get(mods[[j]]$call$x))%in%sub),]
                  test_r[[names(mods)[j]]] <- resp[!(rownames(get(mods[[j]]$call$x))%in%sub), var]
                }else{
                  if ("glm"%in%class(mods[[j]])){
                    train_p[[names(mods)[j]]] <- mods[[j]]$model[sub,2:ncol(mods[[j]]$model),drop=FALSE]
                    train_r[[names(mods)[j]]] <- resp[sub,var]
                    test_p[[names(mods)[j]]] <- mods[[j]]$model[!(rownames(mods[[j]]$model)%in%sub),2:ncol(mods[[j]]$model),drop=FALSE]
                    test_r[[names(mods)[j]]] <- resp[!(rownames(mods[[j]]$model)%in%sub), var]
                  }
                }
              }
            }
            
            
            # check if for all categorical variables condition is met (compare training vs test dat from all models)
            train_comb <- do.call(cbind, train_p)
            test_comb <- do.call(cbind, test_p)
            
            # names of categorical variables
            catVar <- colnames(train_comb)[-which(unlist(lapply(train_comb, is.numeric), use.names = FALSE))]
            levTest <- c()
            for (k in catVar){
              levTest[k] <- all(test_comb[,k]%in%train_comb[,k])
            }
            
            # exit if the condition is met for all categorical variables
            if (all(levTest)) break
          }
          
          # train models on training set
          modsNew <- list()
          for (l in names(mods)){
            if (all(class(mods[[l]]) == "randomForest" )){
              modsNew[[l]] <- randomForest::randomForest(x = train_p[[l]], y = train_r[[l]], ntree = 1000, importance = FALSE)
            }else{
              if("glm"%in%class(mods[[l]])){
                modsNew[[l]] <- glm(train_r[[l]] ~., data = train_p[[l]], family = mods[[l]]$family$family)
              }
            }
          }
          # data frame with observed  for test data set and predicted values from each model
          dat <- data.frame(Observed = test_r[[1]])
          for(m in names(modsNew)){
            if("glm"%in%class(modsNew[[m]])){
              dat[,m] <- predict(modsNew[[m]], newdata = test_p[[m]], type = "response")
            }else{
              dat[,m] <- predict(modsNew[[m]], newdata = test_p[[m]])
            }
          }
          
          # calculate weighted ensemble prediction from individual prediction
          dat[,'ensemblePrediction'] <- rep(10,nrow(dat))
          for (n in 1:nrow(dat)){
            dat[n,'ensemblePrediction'] <- weighted.mean(dat[n,2:(ncol(dat)-1)],weights)
          }
          
          res$cor[i] <- cor(dat$Observed, dat$ensemblePrediction, method = method)
          res$MSE[i] <- mean((dat$Observed - dat$ensemblePrediction)^2)
          
          if(plot == T){
            p <- ggplot(dat, aes(x = Observed, y = ensemblePrediction)) + geom_point() +
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