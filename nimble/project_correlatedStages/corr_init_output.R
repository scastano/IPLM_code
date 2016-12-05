## File called by corr_mcmc.R 

## ####################################
## Generate new directory for output ##
setwd(baseDir)

(currentMcmcDir <- system("find . -type d|grep mcmc", TRUE))
if (length(currentMcmcDir) == 0) {
    system("mkdir mcmc")
}

setwd("mcmc")
(mcmcDir <- paste(format(Sys.time(), "%b%d_%Hh%M.%S"), qsubID, sep="_"))
system( paste("mkdir", mcmcDir, sep=" "))
setwd(mcmcDir)
mcmcDir <- getwd()
nimPrint("MCMC output being saved in directory: ", mcmcDir)


