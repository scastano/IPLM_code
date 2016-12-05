## #########################
## LOAD INTERPOLATION DATA #
setwd(interpDir)

if (CASE== "CLM") {
    ## ######
    ## EGG ## 
    load("interp_cosine_egg_CLM_1000Lines.Rdata")
    if (TempAmpl=="15-25") {
        T <- X.15.25
        listE <- list(a1X.15.25, a2X.15.25, surX.15.25, as.matrix(resVec))
    }## 
    if (TempAmpl=="15-30") {
        T <- X.15.30
        listE <- list(a1X.15.30, a2X.15.30, surX.15.30, as.matrix(resVec))
    }
    ## ########
    ## LARVA ##
    load("interp_cosine_LARVA_CLM_1000Lines.Rdata")
    if (TempAmpl=="15-25") 
        listL <- list(a1X.15.25, a2X.15.25, surX.15.25, as.matrix(resVec))
    ## 
    if (TempAmpl=="15-30") 
        listL <- list(a1X.15.30, a2X.15.30, surX.15.30, as.matrix(resVec))
    ## #######
    ## PUPA ## 
    load("interp_cosine_PUPA_CLM_1000Lines.Rdata")
    if (TempAmpl=="15-25") 
        listP <- list(a1X.15.25, a2X.15.25, surX.15.25, as.matrix(resVec))
    ## 
    if (TempAmpl=="15-30") 
        listP <- list(a1X.15.30, a2X.15.30, surX.15.30, as.matrix(resVec))
    ## #####
    ## GC ##
    load("interp_cosine_gonotrophicCycle_CLM_1000Lines.Rdata")
    if (TempAmpl=="15-25") 
        listGC <- list(a1X.15.25, a2X.15.25, surX.15.25, fecX.15.25, as.matrix(resVec))
    ## 
    if (TempAmpl=="15-30") 
        listGC <- list(a1X.15.30, a2X.15.30, surX.15.30, fecX.15.30, as.matrix(resVec))
}

if (CASE== "IPLM") {
    ## ######
    ## EGG ## 
    load("interp_cosine_egg_IPLM_1000Lines.Rdata")
    if (TempAmpl=="15-25") {
        T <- X.15.25
        listE <- list(a1X.15.25, a2X.15.25, surX.15.25, as.matrix(resVec))
    }## 
    if (TempAmpl=="15-30") {
        T <- X.15.30
        listE <- list(a1X.15.30, a2X.15.30, surX.15.30, as.matrix(resVec))
    }
    ## ########
    ## LARVA ##
    load("interp_cosine_LARVA_IPLM_1000Lines.Rdata")
    if (TempAmpl=="15-25") 
        listL <- list(a1X.15.25, a2X.15.25, surX.15.25, as.matrix(resVec))
    ## 
    if (TempAmpl=="15-30") 
        listL <- list(a1X.15.30, a2X.15.30, surX.15.30, as.matrix(resVec))
    ## #######
    ## PUPA ## 
    load("interp_cosine_PUPA_IPLM_1000Lines.Rdata")
    if (TempAmpl=="15-25") 
        listP <- list(a1X.15.25, a2X.15.25, surX.15.25, as.matrix(resVec))
    ## 
    if (TempAmpl=="15-30") 
        listP <- list(a1X.15.30, a2X.15.30, surX.15.30, as.matrix(resVec))
    ## #####
    ## GC ##
    load("interp_cosine_gonotrophicCycle_IPLM_1000Lines.Rdata")
    if (TempAmpl=="15-25") 
        listGC <- list(a1X.15.25, a2X.15.25, surX.15.25, fecX.15.25, as.matrix(resVec))
    ## 
    if (TempAmpl=="15-30") 
        listGC <- list(a1X.15.30, a2X.15.30, surX.15.30, fecX.15.30, as.matrix(resVec))
}

lT <- length(T) ## 365

## Randomly select nLinesPlot of MCMC output & use it for any stage
indx <-sort(sample(length(resVec), nLinesPlot, rep = (nLinesPlot > length(resVec))))

### ##################
## DEFINE FUNCTION  ##
takeParams2 <- function(indx, listStg, T) { ## X: vector of cosine-profile temperatures over 1 year
    lT <- length(T)
    ## Filter lines
    lList <- length(listStg)
    for (i in 1:lList)
        listStg[[i]] <- listStg[[i]][indx,]
    ##
    Alphas <- cbind(listStg[[1]],listStg[[2]]) ## ncol=2*lT
    lx <- 2*lT
    mu <- t(apply(Alphas, 1, function(x,lx) { x[1:(lx/2)] / (x[1:(lx/2)] + x[(lx/2+1):lx]) }, lx=lx))
    sc <- t(apply(Alphas, 1, function(x,lx) { 1 / (x[1:(lx/2)] + x[(lx/2+1):lx] + 1) }, lx=lx))
    ##
    ## Replace alphas by mu & sc in listStg
    listStg[[1]] <- mu
    listStg[[2]] <- sc
    if(lList > 4) ## case GC
        names(listStg) <- c("mu","sc","sur", "fec", "res") 
    else
        names(listStg) <- c("mu","sc","sur","res") 
    return(listStg)
}


## ############################################
## OBTAIN PARAMETERS FOR EVERY STAGE & CHECK ##
parListE <- takeParams2(indx=indx, listStg=listE, T=T)
class(parListE); length(parListE)
muE  <- parListE[["mu"]]; summary(muE)
scE  <- parListE[["sc"]]; summary(scE)
surE <- parListE[["sur"]]; summary(surE)
resE <- parListE[["res"]]; length(parListE[["res"]]); resE[1:10]


parListL <- takeParams2(indx=indx, listStg=listL, T=T)
class(parListL); length(parListL)
muL  <- parListL[["mu"]]; summary(muL)
scL  <- parListL[["sc"]]; summary(scL)
surL <- parListL[["sur"]]; summary(surL)
resL <- parListL[["res"]]; length(parListL[["res"]]); resL[1:10]


parListP <- takeParams2(indx=indx, listStg=listP, T=T)
class(parListP); length(parListP)
muP  <- parListP[["mu"]]; summary(muP)
scP  <- parListP[["sc"]]; summary(scP)
surP <- parListP[["sur"]]; summary(surP)
resP <- parListP[["res"]]; length(parListP[["res"]]); resP[1:10]


parListGC <- takeParams2(indx=indx, listStg=listGC, T=T)
class(parListGC); length(parListGC)
muGC  <- parListGC[["mu"]]; summary(muGC)
scGC  <- parListGC[["sc"]]; summary(scGC)
surGC <- parListGC[["sur"]]; summary(surGC)
fec   <- parListGC[["fec"]]; summary(fec)
resGC <- parListGC[["res"]]; length(parListGC[["res"]]); resGC[1:10]
