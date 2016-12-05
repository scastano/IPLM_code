## ####################################
## Generate new directory for output ##
setwd(baseDir)

(currentMcmcDir <- system("find . -type d|grep mcmc", TRUE))
if (length(currentMcmcDir) == 0) {
    system("mkdir mcmc")
}

setwd("mcmc")
(mcmcDir <- paste(format(Sys.time(), "%b%d_%Hh%M.%S"), qsubID, resGC, sep="_"))
system( paste("mkdir", mcmcDir, sep=" "))
setwd(mcmcDir)
mcmcDir <- getwd()
nimPrint("MCMC output being saved in directory: ", mcmcDir)


## Copy relevant data,initial cond & parameters info to output directory
system( paste("cp ../../gonoDataConstantsInitial.R ."))
print("gonoDataConstantsInitial.R copied to output directory")

## Initialise output files
mcmcFile       <- "mcmc_samples.csv"
mcmcFile2      <- "mcmc_samples2.csv"
parametersFile <- "mcmc_parameters.txt"




