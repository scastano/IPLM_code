## To correctly initilialise this script, the current directory must be set to a given stage 
## and temps must have been already loaded via (*STG*DataConstantsInitial.R)

setwd(StgDir)
## #############################
## Set directory for plotting ##
## #############################
if (!is.element("figures", dir())) {
    system("mkdir figures")
}
setwd("figures")
getwd()

## ##############
## Set colours ##
rugCol                    <- rgb(1, 0.5, 0.5, 0.8)
nParaCols                 <- 4
nTempCols                 <- length(10:40)
parameterColorsAlpha01    <- rainbow(nParaCols, alpha=0.01)    
parameterColorsAlpha005   <- rainbow(nParaCols, alpha=0.005)    
temperatureColors         <- rainbow(nTempCols)
temperatureColorsAlpha01  <- rainbow(nTempCols, alpha=0.01)    
temperatureColorsAlpha001 <- rainbow(nTempCols, alpha=0.001)
temperatureColorsAlpha005 <- rainbow(nTempCols, alpha=0.05)
iTemp                     <- is.element(10:40, temps)
temperatureColors         <- temperatureColors[iTemp]
temperatureColorsAlpha01  <- temperatureColorsAlpha01[iTemp]
temperatureColorsAlpha001 <- temperatureColorsAlpha001[iTemp] 
temperatureColorsAlpha005 <- temperatureColorsAlpha005[iTemp]
