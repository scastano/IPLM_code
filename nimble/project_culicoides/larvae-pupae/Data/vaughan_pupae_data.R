##############################################################################################
## Vaughan et al's published data of pupal development and mortality at varous temperatures ## 
##############################################################################################

VaughanTemperatures           <- c(20, 23, 27, 30, 35)
VaughanPupaeSurvival          <- c(80.0, 87.7, 94.3, 94.4, 88.9) / 100 ## Published survival proportions
VaughanPreMortalitySampleSize <- 597 ## Total amongst all five groups.
EPreMortalitySampleSize       <- VaughanPreMortalitySampleSize / 5

VaughanPupaeDataHours <- rbind(
    c(temp=20, muOb = 89.7, sigOb = 7.8, ESampleSpre = EPreMortalitySampleSize, ESampleSpost=EPreMortalitySampleSize * VaughanPupaeSurvival[1], Survival=VaughanPupaeSurvival[1]),
    c(temp=23, muOb = 65.5, sigOb = 7.4, ESampleSpre = EPreMortalitySampleSize, ESampleSpost=EPreMortalitySampleSize * VaughanPupaeSurvival[2], Survival=VaughanPupaeSurvival[2]),
    c(temp=27, muOb = 50.6, sigOb = 5.9, ESampleSpre = EPreMortalitySampleSize, ESampleSpost=EPreMortalitySampleSize * VaughanPupaeSurvival[3], Survival=VaughanPupaeSurvival[3]),
    c(temp=30, muOb = 39.1, sigOb = 3.9, ESampleSpre = EPreMortalitySampleSize, ESampleSpost=EPreMortalitySampleSize * VaughanPupaeSurvival[4], Survival=VaughanPupaeSurvival[4]),
    c(temp=35, muOb = 38.8, sigOb = 4.4, ESampleSpre = EPreMortalitySampleSize, ESampleSpost=EPreMortalitySampleSize * VaughanPupaeSurvival[5], Survival=VaughanPupaeSurvival[5])
)

## Convert hours to days 
VaughanPupaeData           <- VaughanPupaeDataHours
VaughanPupaeData[,"muOb"]  <- VaughanPupaeDataHours[,"muOb"] / 24 
VaughanPupaeData[,"sigOb"] <- VaughanPupaeDataHours[,"sigOb"] / 24
VaughanPupaeData           <- as.data.frame(VaughanPupaeData)

## Initial imputation of sample sizes via rejection sampling
## Constraints are (1) the published surviving proportions (2) the total pre-mortality sample size
Precision4ImputingVaughan <- 2 ## 3 ## Poor MCMC mixing when precision = 3. Not an issue when precision = 2. 
imputedDataMatch <- FALSE
iters <- 0
while(!imputedDataMatch) {
    iters <- iters + 1
    (ImpPreMortSS    <- rmultinom(n=1, size=VaughanPreMortalitySampleSize, prob = rep(1/5,5)))
    (ImpPostMortSS   <- round(ImpPreMortSS * VaughanPupaeSurvival))
    imputedDataMatch <- all(round(ImpPostMortSS / ImpPreMortSS, Precision4ImputingVaughan) == round(VaughanPupaeSurvival, Precision4ImputingVaughan))
}
print(paste(iters, "rejection samples required to impute Vaughan's missing sample sizes"))
VaughanPupaeData$ImpPreMortSS  <- ImpPreMortSS
VaughanPupaeData$ImpPostMortSS <- ImpPostMortSS

## Clean up 
rm(iters, ImpPreMortSS, ImpPostMortSS, imputedDataMatch, VaughanPupaeDataHours, EPreMortalitySampleSize, VaughanPreMortalitySampleSize, VaughanTemperatures, VaughanPupaeSurvival)
ls()
VaughanPupaeData
sum(VaughanPupaeData$ImpPreMortSS)

