#######################################################################
## Combined Larvae-> Pupae -> Adult data from Mullens' lab notebooks ##
#######################################################################
options(width=600)

###############
## READ DATA ##
###############
## daysL2A    - number of days since eggs hatched
## repository - experimental 
## male       - number of larvae that became male   adults on daysL2A
## female     - number of larvae that became female adults on daysL2A
## totalAd    - male + female | repostitory
## totalByDay - totalAd summed over repositories
larva17 <- read.table(file="larva-to-adult17.csv", sep=",", header=TRUE)
larva20 <- read.table(file="larva-to-adult20.csv", sep=",", header=TRUE)
larva23 <- read.table(file="larva-to-adult23.csv", sep=",", header=TRUE)       
larva27 <- read.table(file="larva-to-adult27.csv", sep=",", header=TRUE)            
larva30 <- read.table(file="larva-to-adult30.csv", sep=",", header=TRUE)            

########################
## Remove dayL column ##
########################
## DayL was calculated by Mullens as follows
## dayL = observed(daysL2A) - mean(PupalDevTime | Data from another experiment)
larva17 <- larva17[,colnames(larva17)!="dayL"]
larva20 <- larva20[,colnames(larva20)!="dayL"]
larva23 <- larva23[,colnames(larva23)!="dayL"]
larva27 <- larva27[,colnames(larva27)!="dayL"]
larva30 <- larva30[,colnames(larva30)!="dayL"]

(nRepositories17 <- length(unique(larva17$repository)))
(nRepositories20 <- length(unique(larva20$repository)))
(nRepositories23 <- length(unique(larva23$repository)))
(nRepositories27 <- length(unique(larva27$repository)))
(nRepositories30 <- length(unique(larva30$repository)))

sampSizePerRepository <- 80

## Data Checks
head(larva17)
head(larva20)
head(larva23)
head(larva27)
head(larva30)

summary(larva17) ## NAs in male and fremale columns - repressent zeros
summary(larva20) ## NAs in male and fremale columns - repressent zeros
summary(larva23) ## No NAs in male and fremale columns
summary(larva27) ## No NAs in male and fremale columns
summary(larva30) ## No NAs in male and fremale columns

## Replacing the bunch of NA by zeros in the male/female columns
larva17$male   <- replace (larva17$male, is.na(larva17$male), 0)
larva17$female <- replace (larva17$female, is.na(larva17$female), 0)
summary(larva17)
#
larva20$male   <- replace (larva20$male, is.na(larva20$male), 0)
larva20$female <- replace (larva20$female, is.na(larva20$female), 0)
summary(larva20)

## List the data tables
MullensL2A      <- vector("list", 5)
MullensL2A[[1]] <- larva17
MullensL2A[[2]] <- larva20
MullensL2A[[3]] <- larva23
MullensL2A[[4]] <- larva27
MullensL2A[[5]] <- larva30

attributes(MullensL2A[[1]])$Temperature <- 17
attributes(MullensL2A[[2]])$Temperature <- 20
attributes(MullensL2A[[3]])$Temperature <- 23
attributes(MullensL2A[[4]])$Temperature <- 27
attributes(MullensL2A[[5]])$Temperature <- 30

(attributes(MullensL2A[[1]])$SampSize <- nRepositories17 * sampSizePerRepository)
(attributes(MullensL2A[[2]])$SampSize <- nRepositories20 * sampSizePerRepository)
(attributes(MullensL2A[[3]])$SampSize <- nRepositories23 * sampSizePerRepository)
(attributes(MullensL2A[[4]])$SampSize <- nRepositories27 * sampSizePerRepository)
(attributes(MullensL2A[[5]])$SampSize <- nRepositories30 * sampSizePerRepository)

(attributes(MullensL2A[[1]])$PostSampSize <- sum(MullensL2A[[1]]$totalAd))
(attributes(MullensL2A[[2]])$PostSampSize <- sum(MullensL2A[[2]]$totalAd))
(attributes(MullensL2A[[3]])$PostSampSize <- sum(MullensL2A[[3]]$totalAd))
(attributes(MullensL2A[[4]])$PostSampSize <- sum(MullensL2A[[4]]$totalAd))
(attributes(MullensL2A[[5]])$PostSampSize <- sum(MullensL2A[[5]]$totalAd))

## Combine repositories
LP17 <- larva17[!is.na(larva17$totalByDay), c("daysL2A", "totalByDay")]
LP20 <- larva20[!is.na(larva20$totalByDay), c("daysL2A", "totalByDay")]
LP23 <- larva23[!is.na(larva23$totalByDay), c("daysL2A", "totalByDay")]
LP27 <- larva27[!is.na(larva27$totalByDay), c("daysL2A", "totalByDay")]
LP30 <- larva30[!is.na(larva30$totalByDay), c("daysL2A", "totalByDay")]

## Plot cumulative maturation time distributions by temperature and repository
##
if (TRUE) { ## FALSE  ## Toggle this to re-run the plotting
    MX <- max(larva17$daysL2A, larva20$daysL2A, larva23$daysL2A, larva27$daysL2A, larva30$daysL2A)
    for (jj in 1:5) { ## jj = 1
        DATA <- MullensL2A[[jj]]
        DATA$cumTA <- 0
        UN <- unique(DATA$repository)
        (filename <- paste("Mullens_L2A_group", jj,"_T", attributes(DATA)$Temperature, ".png", sep=""))
        png(file=filename)
        for (rp in UN) {   
            index <- DATA$repository == rp
            ss <- DATA[index,] ## a subset of the data for a given repository    
            DATA[index,"cumTA"] <- cumsum(ss$totalAd)
        }
        x <- 0:MX ## max(DATA$daysL2A)
        plot(x, 0*x, ylim=0:1,
             main=paste(attributes(DATA)$Temperature,"degrees"),
             ylab="Proportion Adults", type="n") ## Start a plot for cumulants
        MyCols <- rainbow(length(UN))
        for (rp in UN) { ## rp = 14
            ii <- which(UN==rp)
            index <- DATA$repository == rp
            ss <- DATA[index,] ## a subset of the data for a given repository
            lines(ss$daysL2A, ss$cumTA/tail(ss$cumTA,1), col=MyCols[ii], lwd=3)
        }
        dev.off()
    }
}

##########################
## Non-Developed Vector ##
##########################
(LPsampleSize       <- c(nRepositories17, nRepositories20, nRepositories23, nRepositories27, nRepositories30) * sampSizePerRepository)
(nDeadOrCensorredLP <- LPsampleSize - c(sum(LP17$totalByDay), sum(LP20$totalByDay), sum(LP23$totalByDay), sum(LP27$totalByDay), sum(LP30$totalByDay)))

##############
## Clean Up ##
##############
rm(larva17, larva20, larva23, larva27, larva30, MullensL2A, sampSizePerRepository, LPsampleSize, 
   nRepositories17, nRepositories20, nRepositories23, nRepositories27, nRepositories30)
ls()
## [1] "LP17" "LP20" "LP23" "LP27" "LP30" "LPsampleSize"


## save.image("mullenslarvae-pupae.Rdata")
