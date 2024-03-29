% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VarImp.glm.R
\name{VarImp.glm}
\alias{VarImp.glm}
\title{Variable importance for a glm}
\usage{
VarImp.glm(mod, perm = 100, plot = TRUE, pt = NULL, decreasing = FALSE)
}
\arguments{
\item{mod}{model of class \code{glm}}

\item{perm}{number of permutations}

\item{plot}{logical. create plot for variable importance. (to be
plotted with \code{\link[gridExtra:arrangeGrob]{gridExtra::grid.arrange()}}). Default \code{TRUE}}

\item{pt}{optional. settings for plots. a ggplot \code{theme}}

\item{decreasing}{logical. sort variable importance in decreasing order. Default = \code{FALSE}}
}
\value{
A list with one or two elements
\itemize{
\item \verb{$varImpDat} a data frame with one row for each predictor and two columns (mean and sd variable importance)
\item \verb{$res$varImpPlot} a grob, returned if plot = TRUE
}
}
\description{
Variable importance for a glm
}
\details{
This function calculates variable importance by permuting the predictor variables one at a time and calculating
the \% increase in mean squared error. The function is analogous to the variable importance created by \code{\link[=randomForest]{randomForest()}}.
Plots can be visualized with \code{\link[=grid.arrange]{grid.arrange()}}, e.g. \code{grid.arrange(object$varImpPlot)}
}
\examples{
set.seed(12)
# stepwise glm
full <- glm(MicEnvMod::resp[,2] ~., data = MicEnvMod::pr, family = gaussian)
null <- glm(MicEnvMod::resp[,2] ~1, data = MicEnvMod::pr, family = gaussian)
glmStep <- MASS::stepAIC(full, scope = list(upper = full, lower = null),
direction = 'both')

# Variance importance 
varImp <- VarImp.glm(mod = glmStep, perm = 100, plot = TRUE)
# plot
gridExtra::grid.arrange(varImp$varImpPlot)
}
