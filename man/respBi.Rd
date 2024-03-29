% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/respBi.R
\name{respBi}
\alias{respBi}
\title{Bivariate response plots}
\usage{
respBi(mod, pred, len = 50, pt = NULL, w = 3, h = 3, type = "full")
}
\arguments{
\item{mod}{model for which to generate the response plot.}

\item{pred}{table of predictor variables}

\item{len}{number of data points for which response is predicted. Default \code{50}}

\item{pt}{graphic settings. A ggplot theme can be provided}

\item{w}{width of the plot panel (inch)}

\item{h}{width of the plot panel (inch)}

\item{type}{full model or stepwise selected model (default \code{full})}
}
\value{
list of grobs. heatmaps for response variable as a function of all pairwise combinations of predictors
}
\description{
This function creates bivariate response plots to evaluate how a variable
responds to two predictors at the same time
}
\details{
This function creates response plots by setting all predictor variables to their mean except the two variables of interest and
computes predictions for the response variable.
Plots can be visualized with \code{\link[gridExtra:arrangeGrob]{gridExtra::grid.arrange()}}, e.g. \code{grid.arrange(grobs=object)}
}
\examples{
# for random forest
set.seed(12)
rf <- randomForest::randomForest(x = MicEnvMod::pr, y = resp[,2], ntree = 1000,
 importance = TRUE, nPerm = 100)
respB <- respBi(mod = rf, pred = MicEnvMod::pr)
gridExtra::grid.arrange(grobs = respB, ncol = 3)

# for stepwise glm
full <- glm(MicEnvMod::resp[,2] ~., data = MicEnvMod::pr, family = gaussian)
null <- glm(MicEnvMod::resp[,2] ~1, data = MicEnvMod::pr, family = gaussian)
glmStep <- MASS::stepAIC(full, scope = list(upper = full, lower = null),
direction = 'both')
respM <- respMono(mod = glmStep, pred = MicEnvMod::pr, scale = "fixed")
gridExtra::grid.arrange(grobs = respM, ncol = 7)
}
\references{
Elith, J., Ferrier, S., Huettmann, F., & Leathwick, J. (2005). The evaluation strip:
A new and robust method for plotting predicted responses from species distribution models. \emph{Ecological Modelling},
186(3), 280-289. https://doi.org/10.1016/j.ecolmodel.2004.12.007
}
\seealso{
\code{\link[=respMono]{respMono()}}
}
