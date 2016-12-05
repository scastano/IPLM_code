#################################################################################################
## Supplementary material to "Lefkovitch matrices meet integral projection models: quantifying ##
## the effects of individual heterogeneity in hidden developmental processes"                  ##
#################################################################################################

########################################################
## A list of Nimble functions for Case study I and II ##
########################################################
library(nimble)


###############
## Contents  ##
###############
## roundFun
## roundFun2
## log10
## dLogitUnif
## rLogitUnif
## Registration of dLogitUnif
## nf_muSc2muVar_sub
## nf_muSc2muVar
## nf_muVar2alpha
## nf_muSc2alpha
## nf_alpha2muVar
## nf_muVar2muSc
## nf_alpha2muSc
## nf_getM
## dConditBetaCopula
## rConditBetaCopula
## pConditBetaCopula
## qConditBetaCopula
## Registration of dConditBetaCopula
## nf_getM_given_quality
## nf_setM
## nf_TW1
## nf_TW2
## dJSMD
## rJSMD
## Registration of dJSMD
## dJSMD_condit_quality 
## rJSMD_condit_quality
## Registration of dJSMD_condit_quality
## nf_IPLM
## setM_q



roundFun <- nimbleFunction(
    run = function (x=double(1)) {
        Integer <- round(x)
        returnType(double(1))
        return(Integer)
    }
)


roundFun2 <- nimbleFunction(  
     run = function (x=double(0)) {   
         Integer <- floor(x + 0.5)     
         returnType(integer(0))         
         return(Integer)                 
     }                                    
)  


log10 <- nimbleFunction(
    run = function (x=double(1)) {
        n  <- nimDim(x)[1]
        xx <- numeric(n)
        ##declare(xx, double(1, n))
        L10 <- log(10)
        for (i in 1:n)
            xx[i] <- log(x[i]) / L10
        returnType(double(1))
        return(xx)
    }
)

## ########################################################
## A DISTRIBUTION GIVING THE LOGIT OF A STANDARD UNIFORM ##
## ########################################################
dLogitUnif <- nimbleFunction (
    ## Returns density of x where
    ##                    y ~ Unif(0,1)
    ##                    x = logit(y)
    run = function(x   = double(0),
                   log = integer(0, default = 0)) {
        returnType(double(0))
        y = ilogit(x)
        probX = y * (1 - y) ## Via change of variable
        if (log)
            return(log(probX))
        return(probX)
    }
)

rLogitUnif <- nimbleFunction (
    ## Generates y ~ Unif(0,1)
    ## Returns   x = logit(y)
    run = function(n = integer(0, default=1)
                   ) {
        returnType(double(0))
        if(n != 1) 
            nimPrint("Warning: rLogitUnif only allows n = 1; Using n = 1.\n")
        y <- runif(1,0,1)
        x <- logit(y)
        return(x)
    }
)

registerDistributions(list(dLogitUnif = list(
                               BUGSdist = "dLogitUnif()",
                               discrete = FALSE,
                               types    = c("value=double(0)"), 
                               pqAvail  = FALSE)))


## ############################
## Parameter Transformations ##
## ############################
nf_muSc2muVar_sub <- substitute( nimbleFunction (
    ## Valid for a vector of 2 values
    ## Since .Machine$double.eps does not compile we use substitute / eval here as a work around
    run = function(muSc = double(1), forceBounds = integer(0, default=0)) {
        ## Convert mean and scale to mean and variance of beta distribution
        ## var = scale * mean * (1 - mean)
        mu   <- muSc[1] 
        sc   <- muSc[2]
        if (forceBounds) {
            ## Default is to avoid this behaviour. 
            if (sc <= 0) {
                print("Warning: Ensure scale is in (0, 1). sc = ...")
                print(sc)
                print("Resetting sc to ...")
                sc <- 0 + DOUBLE_EPS
                print(sc)
            } else if (sc >= 1) {
                print("Warning: Ensure scale is in (0, 1). sc = ...")
                print(sc)
                print("Resetting sc to ...")
                sc <- 1 - DOUBLE_EPS
                print(sc)            
            }
        }
        var <- mu * (1-mu) * sc
        ans <- numeric(2)
        ans[1] = mu
        ans[2] = var
        return(ans)
        returnType(double(1))
    }
), list(DOUBLE_EPS = .Machine$double.eps))
## EVALUATE TO OBTAIN NIMBLE FUNCTION
nf_muSc2muVar <- eval(nf_muSc2muVar_sub)



nf_muVar2alpha <- nimbleFunction(
    run = function(muVar = double(1), forceBounds = integer(0, default=0)) {
        ## Valid for a vector of 2 values
        ## Convert mean and variance of beta distribution to canonical parameters
        mu  <- muVar[1]
        var <- muVar[2]
        if (forceBounds) {
            if (var > mu*(1-mu)) {
                print("Warning: var > mu * (1-mu). var & mu*(1-mu) are respectively...")
                print("Resetting so var = mu*(1-mu) = ...")
                var <- mu*(1-mu) ## This maximises the variance, i.e. sets alpha1 = alpha2 = 0 + rounding error
            }
        }
        denom  <- mu * (1 - mu) / var
        alpha1 <- mu * (denom - 1)
        alpha2 <- denom - 1 - alpha1
        ans    <- numeric(2)
        ans[1] <- alpha1
        ans[2] <- alpha2
        return(ans)
        returnType(double(1))
    }
)


nf_muSc2alpha <- nimbleFunction (
    run = function(muSc = double(1)) {
        return(nf_muVar2alpha(nf_muSc2muVar(muSc)))
        returnType(double(1))
    }
)


nf_alpha2muVar <- nimbleFunction (
    run = function(alphas = double(1)) {
        ## Convert canonical parameters of beta distribution to mean and variance
        if ((sum(alphas) < 0))
            stop("Parameters of a beta distribution must be strictly positive.")
        alpha1 <- alphas[1]
        alpha2 <- alphas[2]
        denom  <- alpha1 + alpha2 + 1
        mu     <- alpha1 / (alpha1 + alpha2)
        var    <- mu * (1 - mu) / denom
        muVar  <- numeric(2)
        muVar[1] <- mu
        muVar[2] <- var
        return(muVar)
        returnType(double(1))
    }
)


nf_muVar2muSc <- nimbleFunction ( 
    run = function(muVar = double(1)) {
        ## Convert mean and variance of beta distribution to mu and scale
        mu  <- muVar[1]
        var <- muVar[2]
        if (var > mu*(1-mu)) {
            print("Warning: var > mu * (1-mu). var & mu*(1-mu) are respectively...")
            print("Resetting so var = mu*(1-mu) = ...")
            var <- mu*(1-mu) ## This maximises the variance, i.e. sets alpha1 = alpha2 = 0 + rounding error
        }
        sc      <- var / (mu * (1-mu))        
        muSc    <- numeric(2)
        muSc[1] <- mu
        muSc[2] <- sc
        return(muSc)
        returnType(double(1))
    }
) 



nf_alpha2muSc <- nimbleFunction (
    run = function(alphas = double(1)) {
        return(nf_muVar2muSc(nf_alpha2muVar(alphas)))
        returnType(double(1))
    }
) 


#####################
## Matrix Function ##
#####################
nf_getM <- nimbleFunction(
    ## Returns projection matrix  from eq. 5 (main text)
    run = function(paras = double(1), ## mu, sc, surv,  all in (0,1) 
                   res   = integer(0) ## Resolution of within-stage development
                   ) {
        Dim <- res + 1 ## Matrix dimension
        ## Create a vector to store development level per class
        dev <- numeric(Dim+1)
        for (i in 0:Dim)
            dev[i+1] <- i / Dim      
        ## Initialise the matrix                                   
        M <- matrix(0, nrow = Dim, ncol = Dim)
        for (i in 1:Dim)
            M[i,i] <- 1
        ## Parameterise the first column of M (pdfD) using beta(mu=paras[1], sc=paras[2])
        Alphas <- nf_muSc2alpha(muSc=paras[1:2]) 
        ## Evaluate pdfD - the discretised pdf for development | survival
        pdfD <- numeric(Dim)
        pdfDlower <- pbeta(dev[1], Alphas[1], Alphas[2])
        for(i in 1:Dim) {
            pdfDupper <- pbeta(dev[i+1], Alphas[1], Alphas[2])
            pdfD[i]   <- pdfDupper - pdfDlower
            pdfDlower <- pdfDupper
        } 
        ## Parameterise all other Columns 
        for (j in 1:res) { ## Loop on columns
            M[j:res,j] <- pdfD[1:(Dim-j)]
        }
        ## Columns must sum to 1, so define final row to comply with this constraint
        if (res > 1) {
            for (i in 1:res) {
                if (!(pdfD[Dim-i] == 0 & sum(pdfD[1:(Dim-i)]) > 0)) 
                    M[Dim, i] <- max(0, 1-sum(M[1:res, i]))
            }            
        } else { ## res==1, Dim==2
            M[2,1] <- max(0, 1-M[1,1])
        }
        ## Apply survival - not to last row
        M[1:Dim,1:res] <- paras[3] * M[1:Dim,1:res]
        ## Return matrix 
        return(M)
        returnType(double(2))
    }
)



###########################################################
## BUILDING A USER-DEFINED ConditBetaCopula DISTRIBUTION ##
###########################################################
## STEP1: evaluate density of the ConditBetaCopula distribution 

## A Gaussian copula relates the random variables 'q' (~ unif(0,1)) and 'u2' (~ Beta(alpha1, alpha2) ).
## dConditBetaCopula evaluates the conditional of the marginal-beta | q.
dConditBetaCopula <- nimbleFunction( 
    run = function (x      = double(0), ## Returns the conditional density of x=u2 (WHO'S MARGINAL IS BETA WITH BETAPARAMS)
                    q      = double(0), ## 'u1'. We assume the marginal of u1= q to be U(0,1).
                    rho    = double(0),
                    alphas = double(1),
                    log    = integer(0, default = 0)) { 
        returnType(double(0))
        x1 <- qnorm(q) ## ## set the first variable of the distrib to the qth quantile of a standard normal
        conditionalMean2 <- rho * x1  ## mean of x2 | x1
        conditionalVar2  <- 1-rho*rho ## variance of x2 | x1
        conditionalSD2   <- sqrt(conditionalVar2)
        x2               <- qnorm(pbeta(x, alphas[1], alphas[2])) 
        logDens <- dnorm(x2, conditionalMean2, conditionalSD2, log = TRUE) +    ## f(x2|x1), f =  st. normal 
            dbeta(x, alphas[1], alphas[2], log = TRUE) - dnorm(x2, log = TRUE)  ## Jacobian
        if(log) return(logDens)
        return(exp(logDens))
    }
)


## STEP 2: Random number generation 
## The rNormCopulaBeta function returns draws of the conditional of the marginal-beta | q. 
## We assume q ~ U(0,1), the marginal CDF of our Gaussian copula.
rConditBetaCopula <- nimbleFunction( 
    run = function (n      = integer(0),
                    q      = double(0),
                    rho    = double(0),
                    alphas = double(1)) {
        returnType(double(0))
        if(n != 1) print("rNormCopula only allows n = 1; use n = 1.")
        x1               <- qnorm(q)  ## set x1 to the qth quantile of a standard normal
        conditionalMean2 <- rho * x1  ## mean of x2 given x1
        conditionalVar2  <- 1-rho*rho ## variance of x2 given x1
        conditionalSD2   <- sqrt(conditionalVar2)
        x2               <- rnorm(1, conditionalMean2, conditionalSD2)
        u2               <- pnorm(x2) ## marginal CDF of x2
        beta2            <- qbeta(u2, alphas[1], alphas[2])
        return(beta2)
    }
)

## STEP 3: Calculate cdf of the distribution
pConditBetaCopula <- nimbleFunction( 
    run = function (x      = double(0), ## x (draw from conditBetaCopula) from which I want to know Prob ( u2 < x), space u2
                    q      = double(0),   
                    rho    = double(0),
                    alphas = double(1),
                    log    = integer(0, default = 0)) { 
        returnType(double(0))
        x1 <- qnorm(q) ## set the first variable of the distrib to the qth quantile of a standard normal
        conditionalMean2 <- rho * x1    ## mean of x2 | x1
        conditionalVar2  <- 1 - rho*rho ## variance of x2 | x1
        conditionalSD2   <- sqrt(conditionalVar2)
        u2  <- pbeta(x, alphas[1], alphas[2]) ## turn x back to the 2nd dimension of copula (marginal CDF of beta)
        x2  <- qnorm(u2) ## come back to X1-X2  space (where I apply the bivariate normal)
        cdf <- pnorm(x2, conditionalMean2, conditionalSD2) 
        if(log) return(log(cdf))
        return(cdf)
    }
)

## STEP 4: Calculate the inverse of the cdf of the distribution
qConditBetaCopula <- nimbleFunction( 
    run = function (x      = double(0), ## x (the CDF value from conditBetaCopula) 
                    q      = double(0),   
                    rho    = double(0),
                    alphas = double(1),
                    log    = integer(0, default = 0)) { 
        returnType(double(0))
        x1               <- qnorm(q) ## set the first variable of the distrib to the qth quantile of a standard normal
        conditionalMean2 <- rho * x1 ## mean of x2 | x1
        conditionalVar2  <- 1 - rho*rho ## variance of x2 | x1
        conditionalSD2   <- sqrt(conditionalVar2)
        x2               <- qnorm(x, conditionalMean2, conditionalSD2) 
        u2               <- pnorm(x2)  
        invCdf           <- qbeta(u2, alphas[1], alphas[2]) ## get the u2-value for which the cdf is x
        if(log) return(log(invCdf))
        return(invCdf)
    }
)

## STEP 5 : Register Distribution
registerDistributions(
    list(dConditBetaCopula = list(
             BUGSdist = "dConditBetaCopula(q,rho,alphas)",
             range = c(0, Inf),
             discrete= FALSE,
             types= c("value=double(0)", "q=double(0)","rho=double(0)", "alphas=double(1)"),
             range=c(0,1), 
             pqAvail = FALSE)
))



####################################################
## Matrix Function -individual quality CONSIDERED ##
####################################################
nf_getM_given_quality <- nimbleFunction(
    run = function(paras   = double(1),  ## mu, sc, surv all in (0,1) 
                   res     = integer(0), ## Resolution of within-stage development
                   quality = double(0),  ## Individual quality
                   rho     = double(0)   ## Gaussian correlation parameter - links quality & beta params
                   ) {
        Dim <- res + 1 ## Matrix dimension including transition to following stage
        ## Create a vector to store development level per class
        dev <- numeric(Dim+1)
        for (i in 0:Dim)
            dev[i+1] <- i / Dim 
        ## Initialise the matrix                                   
        M <- matrix(0, nrow = Dim, ncol = Dim)
        for(i in 1:Dim)
            M[i,i] <- 1
        ## Parameterise the first column of M (pdfD) using beta(mu=paras[1], sc=paras[2])
        Alphas <- nf_muSc2alpha(muSc=paras[1:2]) 
        ## Evaluate pdfD - the discretised pdf for development | survival & quality
        pdfD <- numeric(Dim)
        pdfDlower <- pConditBetaCopula(x=dev[1], q=quality, rho=rho, alphas=Alphas) 
        for(i in 1:Dim) {
            pdfDupper <- pConditBetaCopula(x=dev[i+1], q=quality, rho=rho, alphas=Alphas) 
            pdfD[i]   <- pdfDupper - pdfDlower
            pdfDlower <- pdfDupper
        } 
        ## Parameterise all other Columns 
        for (j in 1:res) { ## Loop on columns
            M[j:res,j] <- pdfD[1:(Dim-j)]
        }
        ## Columns must sum to 1, so define final row to comply with this constraint
        if (res > 1) {
            for (i in 1:res) {
                if (!(pdfD[Dim-i] == 0 & sum(pdfD[1:(Dim-i)]) > 0))
                    M[Dim, i] <- max(0, 1-sum(M[1:res, i]))
            }            
        } else { ## res==1, Dim==2
            M[2,1] <- max(0, 1-M[1,1])
        }
        ## Apply survival
        M[1:Dim,1:res] <- paras[3] * M[1:Dim,1:res]
        ## Return matrix given by eq 9 in (25 April version of) paper
        return(M)
        returnType(double(2))
    }
)


nf_setM <- nimbleFunction(
    ## Here we leave all resizing to the nf_IPLM code
    run = function(paras = double(1), ## mu,sc,surv all in (0,1)
                   M     = double(2), 
                   dev   = double(1),
                   pdfD  = double(1)
                   ) {
        Dim <- nimDim(M)[1] 
        Dim1 <- nimDim(dev)[1]
        if (Dim1 != Dim + 1)
            stop("Vector dev not the correct dimension in nf_setM.")
        if (nimDim(pdfD)[1] != Dim)
            stop("Vector pdfD not the correct dimension in nf_setM.")
        res <- Dim-1        ## WITHIN-STAGE resolution.
        ## Fill vector with store development level per class
        for (i in 0:Dim)
            dev[i+1] <- i / Dim  ## dev <- nimVector(0,res)
        ## Initialise the matrix 
        for(i in 1:Dim)
            M[i,i] <- 1
        ## Parameterise the first column of M (pdfD) using beta(mu=paras[1], sc=paras[2])
        Alphas <- nf_muSc2alpha(muSc=paras[1:2]) 
        ## Evaluate pdfD - the discretised pdf for development | surival
        pdfDlower <- pbeta(dev[1], Alphas[1], Alphas[2])
        for(i in 1:Dim) {
            pdfDupper <- pbeta(dev[i+1], Alphas[1], Alphas[2])
            pdfD[i]   <- pdfDupper - pdfDlower
            pdfDlower <- pdfDupper
        } 
        ## Parameterise all other Columns 
        for (j in 1:res) { ## Loop on columns
            M[j:res,j] <- pdfD[1:(Dim-j)]
        }
        ## Columns must sum to 1, so define final row to comply with this constraint
        if (res > 1) {
            for (i in 1:res) {
                if (!(pdfD[Dim-i] == 0  & sum(pdfD[1:(Dim-i)]) > 0)) ## (pdfD[Dim-i] != 0)
                    M[Dim, i] <- max(0, 1-sum(M[1:res, i]))
            }            
        } else { ## res==1, Dim==2
            M[2,1] <- max(0, 1-M[1,1])
        }
        ## Apply survival - not to last row
        M[1:Dim,1:res] <- paras[3] * M[1:Dim,1:res]
        return(M) 
        returnType(double(2))
    }
)
 

#######################
## Travelling Wave 1 ##
#######################
nf_TW1 <- nimbleFunction(
    ## Version of TW used in rJSMD & rJSMD_condit_quality 
    ## dJSMD & dJSMD_condit_quality use nf_TW2
    ## STOPPING RULE: Nstill < thresh | iter == iterMax 
    run = function (M       = double(2),
                    iterMax = integer(0), ## Ten year default maximum
                    thresh  = double(0, default=1E-11) ) { 
        ## Initialise
        Dim       <- nimDim(M)[1]
        res       <- Dim - 1
        ## Nini    <- nimVector(0, Dim)
        Nini      <- numeric(Dim)
        Nini[1]   <- 1        
        N         <- Nini                     ## Initialise population
        End       <- FALSE
        tDevCum   <- numeric(iterMax)         ## Cumulative nb. completing each time step
        tDeadCum  <- tDevCum                  ## Cumulative nb. dead each time step
        iter      <- 0
        while (!End) {
            iter           <- iter + 1
            MN             <- M %*% N          
            N[]            <- MN[,1]           ## for (i in 1:Dim) N[i] <- MN[i,1]
            tDevCum[iter]  <- N[Dim]
            tDeadCum[iter] <- 1 - sum(N)
            Nstill         <- sum(N[1:res])    ## Everybody who's still developing
            End            <- (iter >= iterMax) | (Nstill < thresh)
            if (End) {
                tDevCum  <- tDevCum[1:iter]    ## truncate output vector
                tDeadCum <- tDeadCum[1:iter]   ## truncate output vector                    
            } 
        } 
        ## Probabilities to become fully developed
        pdfDev       <- numeric(iter)
        tDevCumlower <- 0
        for(i in 1:iter) { 
            tDevCumUpper <- tDevCum[i]
            pdfDev[i]    <- tDevCumUpper - tDevCumlower
            tDevCumlower <- tDevCumUpper
        }
        ## Probabilities of dying
        pdfDead       <- numeric(iter)
        tDeadCumlower <- 0
        for(i in 1:iter) {
            tDeadCumUpper <- tDeadCum[i]
            pdfDead[i]    <- tDeadCumUpper - tDeadCumlower
            tDeadCumlower <- tDeadCumUpper
        }
        ## Filter out potential tiny -ves which can occassionaly arise due to rounding error
        for (i in 1:iter) {
            if(pdfDev[i]  < 0) {
                if (pdfDev[i] < -1E-11)  ## Filter rounding error from more serious bugs.
                    nimPrint("Warning: pdfDev[", i, "] = ", pdfDev[i], "\n")
                pdfDev[i]  <- 0
            }
            if(pdfDead[i] < 0) {
                if (pdfDead[i] < -1E-11) ## Filter rounding error from more serious bugs.
                    nimPrint("Warning: pdfDead[", i, "] = ", pdfDead[i], "\n")
                pdfDead[i] <- 0
            }
        }
        ## Combine vectors and return output
        pdfD     <- matrix(0, nrow = iter, ncol = 2)
        pdfD[,1] <- pdfDev
        pdfD[,2] <- pdfDead
        return(pdfD)
        returnType(double(2))
    }
)


#######################
## Travelling Wave 2 ##
#######################
nf_TW2 <- nimbleFunction( ## stop condition : nIter
    run = function (M = double(2), nIter = integer()) {
        Dim       <- nimDim(M)[1]
        res       <- Dim - 1
        ## Initialise population vector
        Nini      <- numeric(Dim)
        Nini[1]   <- 1        
        N         <- Nini ## Initialise population
        End       <- FALSE
        tDevCum   <- numeric(nIter)   ## Cumulative number of individuals completing  stage at each time step
        tDeadCum  <- tDevCum          ## Cumulative number of individuals dead before completing  stage at each time step
        ##
        iter <- 0
        while (!End) {
            iter <- iter + 1
            if (!End) {
                MN             <- M %*% N          
                N[]            <- MN[,1]
                tDevCum[iter]  <- N[Dim]
                tDeadCum[iter] <- 1 - sum(N)
                Condition      <- iter >= nIter
                if (Condition)
                    End <- TRUE
            }      
        }          
        ## Probabilities to become fully developed
        ##pdfDev       <- nimVector(0,nIter)
        pdfDev       <- numeric(nIter)
        tDevCumlower <- 0
        for(i in 1:nIter) { 
            tDevCumUpper <- tDevCum[i]
            pdfDev[i]    <- tDevCumUpper - tDevCumlower
            tDevCumlower <- tDevCumUpper
        }
        ## Probabilities of dying
        ##pdfDead       <- nimVector(0,nIter)
        pdfDead       <- numeric(nIter)
        tDeadCumlower <- 0
        for(i in 1:nIter) {
            tDeadCumUpper <- tDeadCum[i]
            pdfDead[i]    <- tDeadCumUpper - tDeadCumlower
            tDeadCumlower <- tDeadCumUpper
        }
        ## Combine vectors and return output
        pdfD <- matrix(0, nrow = nIter, ncol = 2)
        pdfD[,1] <- pdfDev
        pdfD[,2] <- pdfDead
        return(pdfD)
        returnType(double(2))
    }
)


#####################################
## BUILDING THE dJSMD DISTRIBUTION ##
#####################################
##
## STEP1:  Evaluate density of Joint Sojourn~Mortality Distribution (dJSMD)
## SC: This is fine for indiv with different kernel parameters, othw nf_TW2 output should be an argument here..
dJSMD <- nimbleFunction(
    run = function (x          = double(1),   ## Data vector: EventTime, EventType, Right Censor Indicator
                    paras      = double(1),   ## mu, sc, sur of beta kernel
                    res        = integer(0),  ## Resolution of within stage development
                    CensorTime = integer(0),  ## Right censor starts at CensorTime+1
                    log        = integer(0, default = 0) 
                    ) {
        returnType(double(0))
        if ( min(paras) <= 0 | max(paras) >= 1) { 
            if (log) 
                return(-Inf) ## All parameter domains are (0,1)
            return (0)
        } else {
            Dim <- res + 1
            CT <- CensorTime
            if (CT > 0) {
                nIter <- min(x[1], CT)
                M      <- matrix(0, nrow = Dim, ncol = Dim)
                probTw <- matrix(0, nrow = nIter, ncol = 2)
                M[1:Dim,1:Dim]      <- nf_getM(paras=paras[1:3], res=res)
                probTw[1:nIter,1:2] <- nf_TW2(M=M, nIter=nIter)
                if (x[1] <= CT) {
                    prob <- probTw[x[1],x[2]] ## Not censored
                } else {
                    prob <- max(0, 1-sum(probTw[1:nIter,1:2])) ## Right censor at CensorTime
                }
            } else {
                ## Right censor at 0
                prob <- 1
            }   
            if (log)
                return(log(prob))
            return (prob)
        }
    }
)

## STEP 2 : Random number generation 
rJSMD <- nimbleFunction(     
    run = function (n          = integer(0),  ## n = 1 random variables   
                    paras      = double(1,3), ## mu, sc, sur of beta kernel
                    res        = integer(0),  ## Resolution of within stage development
                    CensorTime = integer(0)   ## Right censor starts at CensorTime+1
                    ) {
        if(n != 1) 
            nimPrint("Warning: rJSMD only allows n = 1; Using n = 1.\n")
        Dim <- res + 1        
        CT  <- CensorTime ## * 1.0 ## Convert to double 
        if (0.0 < CT) {
            M                   <- matrix(0, nrow = Dim, ncol = Dim)
            probTw              <- matrix(0, nrow = CensorTime, ncol = 2)
            pdfD                <- numeric(2*CensorTime + 1)
            M[1:Dim,1:Dim]      <- nf_getM(paras=paras[1:3], res=res)
            probTw[1:CT,1:2]    <- nf_TW2(M=M, nIter=CT) 
            pdfD[1:CT]          <- probTw[,1] ## pDev
            pdfD[(CT+1):(2*CT)] <- probTw[,2] ## pDie
            pdfD[2*CT + 1]      <- max(0, 1-sum(pdfD[1:(2*CT)])) ## pRightCensor
            TimeAndType         <- rcat(1, prob=pdfD)
        } else {
            CT          <- 0 ## A hack to get the following logical arguments working
            TimeAndType <- 1 ## A hack to get the following logical arguments working
        }
        xx <- numeric(3) ## EventTime, EventType, interval
        if (TimeAndType <= CT) {
            ## DEVELOPED
            ## print("Scenario 1")
            xx[1] <- TimeAndType ## EventTime
            xx[2] <- 1           ## EventType = Developed
            xx[3] <- 0 # rinterval(1, t=xx[1], c=CensorTime) ## 0 for t <= c
        } else if (TimeAndType <= 2*CT) {
            ## DIED
            ## print("Scenario 2")
            xx[1] <- TimeAndType - CT ## EventTime
            xx[2] <- 2                ## EventType = Died
            xx[3] <- 0 # rinterval(1, t=xx[1], c=CensorTime) ## 0 for t <= c
        } else {
            ## RIGHT CENSORED
            ## print("Scenario 3")
            xx[1] <- 10^ceiling(log(2+CensorTime)/log(10)) - 1 ## Will return 9, 99, 999, 9999 etc
            xx[2] <- 9
            xx[3] <- 1 # rinterval(1, t=xx[1], c=CensorTime) ## 1 for t > c
        }
        ## Output
        returnType(double(1,3))
        return(xx)
    }
)

## STEP 3 : Register Distribution
registerDistributions(
    list(dJSMD = list(
             BUGSdist = "dJSMD(paras, res, CensorTime)",
             ## range    = c(0, Inf),
             discrete = TRUE,
             types    = c("value=integer(1)","paras=double(1)", "res=integer(0)", "CensorTime=integer(0)"),
             pqAvail  = FALSE))
)



#######################################
## dJSMD_condit_quality DISTRIBUTION ##  
#######################################
dJSMD_condit_quality <- nimbleFunction(
    run = function (x = double(1),           ## Data vector: EventTime, EventType, Right Censor Indicator
                    paras = double(1),       ## mu, sc, sur
                    res        = integer(0), ## Resolution of within stage development
                    CensorTime = integer(0), ## Right censor starts at CensorTime+1
                    quality = double(0),     ## individual quality.
                    rho = double(0),         ## correlation parameter
                    log = integer(0, default = 0)) { 
        returnType(double(0))
        if ( min(paras) <= 0 | max(paras) >= 1) {
            if (log) 
                return(-Inf)
            return (0)
        } else {
            if(res < 1) res   <- 1
            if(res > 100) res <- 100
            Dim <- res + 1
            if (CensorTime > 0) { 
                nIter <- min(x[1], CensorTime)
                M <- matrix(0, nrow = Dim, ncol = Dim)
                probTw <- matrix(0, nrow = nIter, ncol = 2)
                if (quality==0 & rho==0) { ## quality not included, valid when INCLUDE_QUALITY is FALSE in nimbleCode
                    M[1:Dim,1:Dim]      <- nf_getM(paras=paras[1:3], res=res)
                    probTw[1:nIter,1:2] <- nf_TW2(M=M, nIter=nIter)
                    ## Not censored;
                    if (x[1] <= CensorTime) 
                        prob <- probTw[x[1],x[2]]
                    ## else : Right censor at CensorTime
                    prob <- max(0, 1-sum(probTw[1:nIter,1:2]))
                } else {   ## Quality included:
                    M[1:Dim,1:Dim]      <- nf_getM_given_quality(paras=paras[1:3], res=res, quality=quality, rho=rho)
                    probTw[1:nIter,1:2] <- nf_TW2(M=M, nIter=nIter)
                    ## Not censored;
                    if (x[1] <= CensorTime) {
                        prob <- probTw[x[1],x[2]]
                    } else {  ## right censor at CensorTime
                        prob <- max(0, 1-sum(probTw[1:nIter,1:2]))
                    }
                }
            } else { ##CensorTime = 0
                prob <- 1 ## individual not followed anymore (either died or censored before)
            } 
            if (log)
                return(log(prob))
            return (prob)
        }
    }
)

## random generator
rJSMD_condit_quality <- nimbleFunction(     
    run = function (n = integer(0),          ## Nb. of required samples
                    paras=double(1),         ## mu, sc, sur
                    res=integer(0),
                    CensorTime = integer(0), ## Right censor starts at CensorTime+1
                    quality = double(0),     ## individual quality.
                    rho = double(0) ){       ## correlation parameter
        if(n != 1) print("Warning: rJSMD only allows n = 1; using n = 1.")
        if(res < 1) res   <- 1
        if(res > 100) res <- 100
        Dim <- res + 1
        if (0.0 < CensorTime) {
            M      <- matrix(0, nrow = Dim, ncol = Dim)
            probTw <- matrix(0, nrow = CensorTime, ncol = 2)
            pdfD   <- numeric(2*CensorTime + 1)
            if (quality==0 & rho==0) { 
                M[1:Dim,1:Dim] <- nf_getM(paras=paras[1:3], res=res)## TRUE when INCLUDE_QUALITY is FALSE in nimbleCode
            } else { ## quality INCLUDED:
                M[1:Dim,1:Dim] <- nf_getM_given_quality(paras=paras[1:3], res=res, quality=quality, rho=rho)
            }
            probTw[1:CensorTime,1:2]    <- nf_TW2(M=M, nIter=CensorTime) 
            pdfD[1:CensorTime]          <- probTw[,1]         ## pDev
            pdfD[(CensorTime+1):(2*CensorTime)] <- probTw[,2] ## pDie
            pdfD[2*CensorTime + 1]      <- max(0, 1-sum(pdfD[1:(2*CensorTime)])) ## pRightCensor
            TimeAndType                 <- rcat(1, prob=pdfD)
        } else { ## CensorTime==0
            CensorTime   <- 0 
            TimeAndType  <- 1 
        }
        xx <- numeric(3) ## EventTime, EventType, interval
        if (TimeAndType <= CensorTime) {
            ## DEVELOPED
            xx[1] <- TimeAndType                         ## EventTime
            xx[2] <- 1                                   ## EventType = Developed
            xx[3] <- 0 
        } else if (TimeAndType <= 2*CensorTime) {
            ## DIED
            xx[1] <- TimeAndType - CensorTime                    ## EventTime
            xx[2] <- 2                                   ## EventType = Died
            xx[3] <- 0 
        } else {
            ## RIGHT CENSORED
            xx[1] <- 10^ceiling(log(2+CensorTime)/log(10)) - 1 ## Will return 9, 99, 999, 9999 etc
            xx[2] <- 9
            xx[3] <- 1 
        }
        ## Output
        returnType(double(1,3))
        return(xx)
    }
)

## STEP 3 : Register Distribution
registerDistributions(
    list(dJSMD_condit_quality = list(
             BUGSdist = "dJSMD_condit_quality(paras, res, CensorTime, quality, rho)",
             range = c(0, Inf),
             discrete= TRUE,
             types= c("value=integer(1)", "paras=double(1)", "res=integer(0)", "CensorTime=integer(0)",
                      "quality=double(0)", "rho=double(0)"),
             pqAvail = FALSE))
)


##############################################
## Build time-distributed Lefkovitch matrix ##
##############################################
nf_IPLM <- nimbleFunction(
    run = function(paras  = double(2), ## mu, sc, surv \in (0,1) X nStage
                   res    = double(1), ## DON'T SET THIS TO interger
                   femfec = double(0), ## Expected female-only fecundity 
                   gCycle = integer(0, default=0) ## 1==Include Gonotrophic cycle
                   ) {
        rowsParas <- nimDim(paras)[1]
        colsParas <- nimDim(paras)[2]
        nStage    <- nimDim(res)[1]
        if (rowsParas != nStage) { nimPrint("Warning: Dimension problems in nf_IPLM\n") }
        DimTotal <- sum(res) 
        ## Initialise the matrices and vectors     
        M <- matrix(0,nrow = DimTotal, ncol = DimTotal) 
        ## Loop on stages
        baseIndex <- 0
        for (i in 1:nStage) { ## i=1
            Resi  <- res[i]
            Dimi   <- Resi + 1
            Dimi1  <- Resi + 2
            PiTi01 <- matrix(0, nrow = Dimi, ncol = Dimi) ##submatrices P - T - 0 - 1 (look at the paper)
            colIndices <- numeric(Resi)
            Ti1 <- numeric(Resi)
            dev    <- numeric(Dimi1)
            pdfD   <- numeric(Dimi)
            PiTi01 <- nf_setM(paras = paras[i,1:3], M = 0*PiTi01, 
                              dev   = 0*dev, pdfD = 0*pdfD)
            for (j in 1:Resi) {
                colIndices[j] <- baseIndex + j
                Ti1[j]        <- PiTi01[Dimi, j]
            }
            M[colIndices[1]:colIndices[Resi],
              colIndices[1]:colIndices[Resi]] <- PiTi01[1:Resi,1:Resi]
            if (i < nStage) {
                M[colIndices[Resi]+1, colIndices[1]:colIndices[Resi]] <- Ti1
            } else {
                M[1,colIndices[1]:colIndices[Resi]] <- Ti1 * femfec
                if (gCycle == 1) {
                    M[colIndices[1], colIndices[1]:colIndices[Resi]] <- Ti1 +
                                      M[colIndices[1], colIndices[1]:colIndices[Resi]]
                }
            }            
            baseIndex <- sum(res[1:i])
        }
        ## Output
        returnType(double(2))
        return(M) 
    }
)


setM_q <- nimbleFunction( ## quality included
    ## nf_getM doesn't like being resized once compiled
    ## Here we leave all resizing to the nf_IPLM code
    run = function(paras   = double(1), ## mu,sc,surv all in (0,1)
                   quality = double(0),
                   rho     = double(0),
                   M       = double(2), ## Matrix dimension when including transition to following stage
                   dev     = double(1),
                   pdfD    = double(1)
                   ) {
        Dim  <- nimDim(M)[1] 
        Dim1 <- nimDim(dev)[1]
        if (Dim1 != Dim + 1)
            stop("Vector dev not the correct dimension in setM_q.")
        if (nimDim(pdfD)[1] != Dim)
            stop("Vector pdfD not the correct dimension in setM_q.")
        res <- Dim-1 ## WITHIN-STAGE resolution.
        ## Fill vector with store development level per class
        for (i in 0:Dim)
            dev[i+1] <- i / Dim
        ## Initialise the matrix 
        for(i in 1:Dim)
            M[i,i] <- 1
        ## Parameterise the first column of M (pdfD) using beta(mu=paras[1], sc=paras[2])
        Alphas     <- nf_muSc2alpha(muSc=paras[1:2]) 
        ## Evaluate pdfD - the discretised pdf for development | survival
        pdfDlower     <- pConditBetaCopula(x=dev[1], q=quality, rho=rho, alphas=Alphas)   
        for(i in 1:Dim) {
            pdfDupper <- pConditBetaCopula(x=dev[i+1], q=quality, rho=rho, alphas=Alphas)
            pdfD[i]   <- pdfDupper - pdfDlower
            pdfDlower <- pdfDupper
        } 
        ## Parameterise all other Columns 
        for (j in 1:res) { ## Loop on columns
            M[j:res,j] <- pdfD[1:(Dim-j)]
        }
        ## Columns must sum to 1, so define final row to comply with this constraint
        if (res > 1) {
            for (i in 1:res) {
                if (pdfD[Dim-i] != 0)
                    M[Dim, i] <- max(0, 1-sum(M[1:res, i]))
            }            
        } else { ## res==1, Dim==2
            M[2,1] <- max(0, 1-M[1,1])
        }
        ## Apply survival - not to last row
        M[1:Dim,1:res] <- paras[3] * M[1:Dim,1:res]
        return(M) 
        returnType(double(2))
    }
)


