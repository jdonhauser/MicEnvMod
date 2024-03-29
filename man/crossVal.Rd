% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossVal.R
\name{crossVal}
\alias{crossVal}
\title{Cross validation by repeated split sampling}
\usage{
crossVal(
  pred,
  resp,
  var,
  rounds = 200,
  frac = 0.7,
  plot = FALSE,
  method = "pearson",
  type = "rf",
  family = "gaussian",
  ignoreError = FALSE
)
}
\arguments{
\item{pred}{data frame or matrix of predictor variables}

\item{resp}{data frame with response variables, rownames have to correspond
to \code{pred}}

\item{var}{response variable (column name or index of \code{resp})}

\item{rounds}{number of split sampling runs. Default \code{200}}

\item{frac}{fraction of samples used for training set. Default \code{0.7}}

\item{plot}{logical. plot observed vs predicted (as list of grobs to be
plotted with grid.arrange()). Default \code{FALSE}}

\item{method}{method for correlation between observed and predicted, as in
\code{\link[=cor]{cor()}}. Default \code{pearson}}

\item{type}{model type; "rf" for random forest or "glm"}

\item{family}{family for glm as in \code{\link[=glm]{glm()}}. Default \code{"gaussian"}}

\item{ignoreError}{logical. If \code{TRUE} an iteration with error is skipped and a new subsample is taken; proceeded until desired number of repeats is reached;
useful e.g for non-gaussian glms that do not converge for the training set. Default \code{FALSE}}
}
\value{
A list with three or four elements
\itemize{
\item \verb{$cor} correlations observed vs predicted
\item \verb{$MSE} Mean squared errors observed vs predicted:
\item \verb{$summary} data frame with mean and sd of all cross validation runs
\item \verb{$plots} plots (list of grobs)
}
}
\description{
This function evaluates model performance by repeated split sampling based on correlations
between observed and predicted values.
}
\details{
This function splits the data in a test and training set by randomly drawing samples without replacement. The function builds the model with the training set,
predicts the response variable for the test test. Measures of accuracy include correlation between observed and predicted values as well as mean squared errors. The procedure is repeated n times
and results are averaged. Measures of accuracy can be used as weights in an ensemble model.
Plots of observed vs predicted values can be visualized with \code{\link[gridExtra:arrangeGrob]{gridExtra::grid.arrange()}}, e.g. \code{grid.arrange(grobs=object$plots)}
}
\examples{
# for random forest
set.seed(12)
rf <- randomForest::randomForest(x = MicEnvMod::pr, y = resp[,2],
ntree = 1000, importance = TRUE, nPerm = 100)
cv <- crossVal(pred = MicEnvMod::pr, resp = MicEnvMod::resp,
var = "optimum_temperature", plot = TRUE, type = "rf")
# summary
cv$summary
# plot
gridExtra::grid.arrange(grobs = cv$plots, ncol = 15)

# for stepwise glm
full <- glm(MicEnvMod::resp[,2] ~., data = MicEnvMod::pr, family = gaussian)
null <- glm(MicEnvMod::resp[,2] ~1, data = MicEnvMod::pr, family = gaussian)
glmStep <- MASS::stepAIC(full, scope = list(upper = full, lower = null), direction = 'both')
set.seed(12)
cv <- crossVal(pred = glmStep$model[,2:ncol(glmStep$model)], resp = MicEnvMod::resp,
var = "optimum_temperature", plot = TRUE, type = "glm")
# summary
cv$summary
# plot
gridExtra::grid.arrange(grobs = cv$plots, ncol = 15)
}
\references{
Thuiller, W., Lafourcade, B., Engler, R., & Araujo, M. B. (2009).
BIOMOD - A platform for ensemble forecasting of species distributions.
\emph{Ecography}, 32(3), 369-373. https://doi.org/10.1111/j.1600-0587.2008.05742.x
}
\seealso{
\itemize{
\item \code{\link[=crossVal.step]{crossVal.step()}} for cross validation with stepwise selection
\item \code{\link[=crossVal.ensemble]{crossVal.ensemble()}} for cross validation of a weighted ensemble model
}
}
