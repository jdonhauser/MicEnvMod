# MicEnvMod: An R package for modelling microbial community - environment relationships

## Installation

```
devtools::install_github("jdonhauser/MicEnvMod")
```
## Functions
### varImp.glm 
Variable importance for a glm by permuting predictor variables
### varImpPlotgg
Plot variable importance from `randomForest` as barplot with ggplot
### respMono
Plot response a dependent variable to a particular predictor. Calculated by setting all predictor variables except one to a fixed value.
### respBi
Bivariate response plot, same as respMono, but with two predictor variables.
### crossVal
model cross validation by repeated split sampling. Currently models of class `randomForest` and `glm` are supported.
### crossVal.step
glm cross validation by repeated split sampling with stepwise selection. Assesses stability of predictor variable selection.
### crossVal.ensemble
Cross validation of weighte ensemble models from different model types by repeated split sampling

## Reference
Donhauser, J., Doménech-Pascual, A., Han, X., Jordaan, K., Ramond, J.-B., Frossard, A., Romaní, A.M., Priemé, A., 2024. Modelling soil prokaryotic traits across environments with the trait sequence database *ampliconTraits* and the R package *MicEnvMod*. Ecological Informatics 83, 102817. doi:10.1016/j.ecoinf.2024.102817
