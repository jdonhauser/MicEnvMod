% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/respMono.R
\name{respMono}
\alias{respMono}
\title{Monovariate response plots}
\usage{
respMono(
  mod,
  pred,
  len = 200,
  pt = NULL,
  w = 2,
  h = 2,
  type = "full",
  scale = "free",
  cols = NULL
)
}
\arguments{
\item{mod}{model for which to generate the response plot.}

\item{pred}{table of predictor variables}

\item{len}{number of data points for which response is predicted. Default \code{200}}

\item{pt}{optional graphic settings. A ggplot theme can be provided}

\item{w}{width of the plot panel (inch)}

\item{h}{width of the plot panel (inch)}

\item{type}{full model or stepwise selected model (default \code{full})}

\item{scale}{\code{"free"} or \code{"fixed"}. if set to \code{"free"} (default), the range of the y axis
(response variable) varies across plots. if set to \code{"fixed"} the same scale is
used for all plots}

\item{cols}{optional named color vector for bars representing levels of a
categorical variable}
}
\value{
list of grobs, scatter plots for continuous variables and barplots for categorical variables as function of each predictor variable
}
\description{
This function creates monovariate response plots to evaluate how a variable
responds to a particular predictor
}
\details{
This function creates response plots by setting all predictor variables to their mean (for continuous variable) or the first level, except the variable of interest and
computes predictions for the response variable.
Plots can be visualized with \code{\link[gridExtra:arrangeGrob]{gridExtra::grid.arrange()}}, e.g. \code{grid.arrange(grobs=object)}
}
\examples{
# for random forest
set.seed(12)
rf <- randomForest::randomForest(x = MicEnvMod::pr, y = resp[,2], ntree = 1000,
 importance = TRUE, nPerm = 100)
respM <- respMono(mod = rf, pred = MicEnvMod::pr, scale = "fixed")
gridExtra::grid.arrange(grobs = respM, ncol = 3)

# for stepwise glm
full <- glm(MicEnvMod::resp[,2] ~., data = MicEnvMod::pr, family = gaussian)
null <- glm(MicEnvMod::resp[,2] ~1, data = MicEnvMod::pr, family = gaussian)
glmStep <- MASS::stepAIC(full, scope = list(upper = full, lower = null),
direction = 'both')
respM <- respMono(mod = glmStep, pred = MicEnvMod::pr, scale = "fixed")
gridExtra::grid.arrange(grobs = respM, ncol = 3)
}
\references{
Elith, J., Ferrier, S., Huettmann, F., & Leathwick, J. (2005). The evaluation strip:
A new and robust method for plotting predicted responses from species distribution models. \emph{Ecological Modelling},
186(3), 280-289. https://doi.org/10.1016/j.ecolmodel.2004.12.007
}
