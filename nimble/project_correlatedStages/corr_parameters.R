## File called by corr_mcmc.R

## #############################
## PARAMETERS for corr_mcmc.R ##
## #############################

setwd("~/IPLM_paper/nimble")
library(nimble)
source("FUNCTIONS_NIMBLE.R")


## CONSTANTS
N          <- 50  ## Sample size
S          <- 2   ## Number of stages
ResMax     <- 100 ## Maximum resolution 
CensorTime <- 21
(Constants <- list(N=N, S=S, ResMax = ResMax, CensorTime=CensorTime))


## PROPOSE INITIAL VALUES
Al <-  matrix(rep(0.5,6),2,3) ## an arbitrary matrix of alphas that violates the constraint
ll <- 0
while(any(Al<=1)) { ## 
    ll <- ll+1
    print(ll)
    kParas1 <- runif(3) ## kernel parameters (mu, sc, sur) for Stage 1
    kParas2 <- runif(3) ## kernel parameters (mu, sc, sur) for Stage 2
    (paras  <- rbind(kParas1,kParas2)) 
    ## CHECKING mu, sc respect the constraint on alpha1, alpha2
    Bpar <- paras[,-3]
    Al   <- apply(Bpar,1,function(x) nf_muVar2alpha(nf_muSc2muVar(x)))
}

## Visual check
if (FALSE) { ## TRUE
    par(mfrow=c(2,1))
    curve(dbeta(x, Al[1,1], Al[1,2]),0,1,n=1001)
    curve(dbeta(x, Al[2,1], Al[2,2]),0,1,n=1001)
}   

(rho     <-runif(1))
(qual    <- sort(runif(N)))
(res     <- sample(1:100,2, replace=TRUE))
(resCont <- res / ResMax)
(res     <- ceiling(ResMax * resCont))

Inits <- list(paras           = paras,
                rho           = rho,
                qual          = qual,
                logit_resCont = logit(resCont),
                resCont       = resCont,
                res           = res,
                logit_paras   = logit(paras),
                logit_rho     = logit(rho),
                logit_qual    = logit(qual))



Data <- list(constraintData=c(1,1), y=array(dim=c(N,3,S)))
