## #########################
## LOAD INTERPOLATION DATA #
setwd(interpDir)


if (CASE== "CLM") {
    ## EGG
    load("interp_egg_CLM_1000Lines.Rdata") 
    listE <- list(a1X, a2X, surX, as.matrix(resVec))
    ## GC
    load("interp_gonotrophicCycle_CLM_1000Lines.Rdata")
    listGC <- list(a1X, a2X, surX, fecX, as.matrix(resVec))
    ## LARVA
    load("interp_LARVA_CLM_1000Lines.Rdata")
    listL <- list(a1X, a2X, surX, as.matrix(resVec))
    ## PUPA
    load("interp_PUPA_CLM_1000Lines.Rdata")
    listP <- list(a1X, a2X, surX, as.matrix(resVec))
}


if (CASE== "IPLM") {
    ## EGG
    load("interp_egg_IPLM_1000Lines.Rdata") 
    listE <- list(a1X, a2X, surX, as.matrix(resVec))
    ## GC
    load("interp_gonotrophicCycle_IPLM_1000Lines.Rdata")
    listGC <- list(a1X, a2X, surX, fecX, as.matrix(resVec)) 
    ## LARVA
    load("interp_LARVA_IPLM_1000Lines.Rdata")
    listL <- list(a1X, a2X, surX, as.matrix(resVec))
    ## PUPA
    load("interp_PUPA_IPLM_1000Lines.Rdata")
    listP <- list(a1X, a2X, surX, as.matrix(resVec))
}


if (CASE== "MAP") { 
    ## EGG
    load("interp_egg_MAP_1000Lines.Rdata") 
    listE <- list(a1X, a2X, surX, as.matrix(resVec))
    ## GC
    load("interp_gonotrophicCycle_MAP_1000Lines.Rdata")
    listGC <- list(a1X, a2X, surX, fecX, as.matrix(resVec)) 
    ## LARVA
    load("interp_LARVA_MAP_1000Lines.Rdata")
    listL <- list(a1X, a2X, surX, as.matrix(resVec))
    ## PUPA
    load("interp_PUPA_MAP_1000Lines.Rdata")
    listP <- list(a1X, a2X, surX, as.matrix(resVec))
}


## ############################################
## Set temperature set to study trajectories ## 
if (PLOT == "AnnualTraj" | PLOT == "SteadyState" | PLOT == "Damping") 
    T <- seq(10, 35, by=5)
if (PLOT == "AnnualTraj" | PLOT == "Damping") 
    T <- seq(10, 35, by=5)
if (PLOT == "PopGrowthRate") 
    T <- seq(10, 40, by=1)

(lT <- length(T))

## Randomly select nLinesPlot of MCMC output & use it for any stage
indx <-sort(sample(length(resVec), nLinesPlot, rep = (nLinesPlot > length(resVec)) ))
## indx <-sort(sample(nrow(survX), nLinesPlot, rep = (nLinesPlot > nrow(survX)) ))

### ##################
## DEFINE FUNCTION  ##
takeParams <- function(indx, listStg, X, T) { ## X: vector of temperatures of interpolated data (0 to 50 by 0.2)
    lT <- length(T); lX <- length(X)
    ## Filter lines
    lList <- length(listStg)
    for (i in 1:lList)
        listStg[[i]] <- listStg[[i]][indx,]
    ## Filter columns for all but the last element of list
    newl <- list()
    for (i in 1:(lList-1)) {
        tmp <- matrix(0,length(indx),lT)
        for (j in 1:lT) {
            (jj <- which(X == T[j])) 
            tmp[,j] <- listStg[[i]][,jj]
        }
        newl[[i]] <- tmp 
    }
    ## Use first two elements of newl (the alphas) to calculate mean and sc:
    Alphas <- cbind(newl[[1]],newl[[2]])
    lx <- 2*lT
    mu <- t(apply(Alphas, 1, function(x,lx) { x[1:(lx/2)] / (x[1:(lx/2)] + x[(lx/2+1):lx]) }, lx=lx))
    sc <- t(apply(Alphas, 1, function(x,lx) { 1 / (x[1:(lx/2)] + x[(lx/2+1):lx] + 1) }, lx=lx))
    ##
    if(lList>4) ## case GC
        flist <- list(mu=mu, sc=sc, sur=newl[[lList-2]], fec=newl[[lList-1]], res=listStg[[lList]])
    else 
        flist <- list(mu=mu, sc=sc, sur=newl[[lList-1]], res=listStg[[lList]])
    ##
    return(flist)
}

## ############################################
## OBTAIN PARAMETERS FOR EVERY STAGE & CHECK ##
parListE <- takeParams(indx=indx, listStg=listE, X=X, T=T)
class(parListE); length(parListE)
muE  <- parListE[["mu"]]; summary(muE)
scE  <- parListE[["sc"]]; summary(scE)
surE <- parListE[["sur"]]; summary(surE)
resE <- parListE[["res"]]; length(parListE[["res"]]); resE[1:10]

parListL <- takeParams(indx=indx, listStg=listL, X=X, T=T)
class(parListL); length(parListL)
muL  <- parListL[["mu"]]; summary(muL)
scL  <- parListL[["sc"]]; summary(scL)
surL <- parListL[["sur"]]; summary(surL)
resL <- parListL[["res"]]; length(parListL[["res"]]); resL[1:10]

parListP <- takeParams(indx=indx, listStg=listP, X=X, T=T)
class(parListP); length(parListP)
muP  <- parListP[["mu"]]; summary(muP)
scP  <- parListP[["sc"]]; summary(scP)
surP <- parListP[["sur"]]; summary(surP)
resP <- parListP[["res"]]; length(parListP[["res"]]); resP[1:10]

parListGC <- takeParams(indx=indx, listStg=listGC, X=X, T=T)
class(parListGC); length(parListGC)
muGC  <- parListGC[["mu"]]; summary(muGC)
scGC  <- parListGC[["sc"]]; summary(scGC)
surGC <- parListGC[["sur"]]; summary(surGC)
fec   <- parListGC[["fec"]]; summary(fec)
resGC <- parListGC[["res"]]; length(parListGC[["res"]]); resGC[1:10]
