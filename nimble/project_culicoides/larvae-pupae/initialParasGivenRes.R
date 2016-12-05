##########################################################################################
## Here we set initial parameter estimates given resL and resP so that LogProb != -Inf  ##
##########################################################################################

oriResL        <- LPModel$resL 
oriResP        <- LPModel$resP 
oriLogitParasL <- LPModel$logit_parasL 
oriLogitParasP <- LPModel$logit_parasP
rL             <- LPModel$resL 
rP             <- LPModel$resP 


## These loops were used to develop a scheme to modify parameter starting values, given res, to avoid -Inf LogProb
LPModel$logit_parasL <- oriLogitParasL
LPModel$logit_parasP <- oriLogitParasP
(LPModel$resL <- rL)
(LPModel$resP <- rP)
rL; rP
LPModel$logit_parasP[1,3] <- oriLogitParasP[1,3] + 2
LPModel$logit_parasL[,1]  <- oriLogitParasL[,1]  + (LPModel$resL<=2) + (LPModel$resL<2)
if (rL>=11) {
    LPModel$logit_parasL[1,1] <- oriLogitParasL[1,1] - 0.08 * (rL>=11&rP>2) - 0.08 * (rL>=12&rP>2) -
        0.2 * (rL>=13&rP>2) - 0.1*(rL>=16&rP>2) - 0.1*(rL>=18&rP>2) - 0.2*(rL>=19&rP>2) - 0.1*(rL>=21) - 0.1*(rL>=27) - 0.1*(rL>=34) - 0.1*(rL>=43) 
    LPModel$logit_parasL[2,1] <- oriLogitParasL[2,1] - 0.08 * (rL>=11&rP>2) - 0.03 * (rL>=12&rP>2) - 0.1 * (rL>=13&rP>2) - 
        0.1*(rL>=14&rP>2) - 0.1*(rL>=17&rP>2) - 0.1*(rL>=19&rP>2) - 0.1*(rL>=21) - 0.1*(rL>=27) - 0.1*(rL>=34) - 0.1*(rL>=43) 
    LPModel$logit_parasL[3,1] <- oriLogitParasL[3,1] - 0.1 * (rL>=14) * min(3,rP) - 0.1*(rL>=21) - 0.1*(rL>=27) - 0.1*(rL>=34)
}
if (rL>=14 & rP>=2) {
    (LPModel$logit_parasL[,1] <- c(-2.8, -2.75, -2.7, -2.32, -2.33))
    (LPModel$logit_parasL[,2] <- LPModel$logit_parasL[,2] - 0.5)
}
if (rL>=20 & rP>2) {
    LPModel$logit_parasL[,1] <- LPModel$logit_parasL[,1] - 0.05
}
if (rL>=21) {
    LPModel$logit_parasL[c(2,4),1] <- LPModel$logit_parasL[c(2,4),1] - c(0.06, 0.03)
    if (rP >= 2) LPModel$logit_parasL[1,1] <- LPModel$logit_parasL[1,1] -0.02
    if (rP >= 5) LPModel$logit_parasL[1,1] <- LPModel$logit_parasL[1,1] - 0.05 
}
if (rL>=22) { if (rP >= 4) LPModel$logit_parasL[1,1] <- LPModel$logit_parasL[1,1] - 0.05 }
if (rL>=23) { if (rP >= 4) LPModel$logit_parasL[1,1] <- LPModel$logit_parasL[1,1] - 0.05 }
if (rL>=25) { if (rP >= 4) LPModel$logit_parasL[1,1] <- LPModel$logit_parasL[1,1] - 0.05 }
if (rL>=27) { if (rP >= 3) LPModel$logit_parasL[1:2,1] <- LPModel$logit_parasL[1:2,1] - 0.05 }
if (rL>=29) { if (rP >= 4) LPModel$logit_parasL[1,1] <- LPModel$logit_parasL[1,1] - 0.05 }
if (rL>=30) { if (rP >= 3) LPModel$logit_parasL[1:2,1] <- LPModel$logit_parasL[1:2,1] - c(0.5, 0.05) }
if (rL>=33) { if (rP >= 3) LPModel$logit_parasL[2,2] <-LPModel$logit_parasL[2,2] + 0.5 }
if (rL>=41) {
    if (rP >= 3) {
        LPModel$logit_parasL[1:2,1] <-LPModel$logit_parasL[1:2,1] - c(0.2, 0.1)
        LPModel$logit_parasL[2,2] <-LPModel$logit_parasL[2,2] + 0.1
    }
} 
rL; rP 
LPModel$logit_parasL
dUnimodal(LPModel$logit_parasL[1:nTempsL, 1], log=1) 
(TEST_LogProb <- LPModel$calculate())
LPModel$logProb_constraintDataL
LPModel$logProb_constraintDataP
if (!is.finite(TEST_LogProb )) {
    nimPrint(rL, " ", rP)
    stop()
} 
nimPrint(paste(rL, rP, TEST_LogProb ))


##########
## TEST ##
(TEST_LogProb  <- LPModel$calculate()); LPModel$logProb_constraintDataL
if (!is.finite(TEST_LogProb ))
    stop("LogProb = -Inf, you must reset some of the starting values.")
nimPrint(paste(rL, rP, TEST_LogProb))

rm(oriResL, oriResP, oriLogitParasL, oriLogitParasP, rL, rP)

