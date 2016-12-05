################################################################################################################################################
## Data from Mullen's culicoides development studies - combines data as published with suplimentary information taken from hand-written notes ## 
## SampleSizePre and SampleSizePost were obtained from the Mullens' paper notes                                                               ##
################################################################################################################################################

MullensPupaeDataHours <- rbind(c(temp=17, muOb = 119.4, sigOb = 10.1, SampleSizePre =  81, SampleSizePost =  78),
                               c(temp=20, muOb =  85.0, sigOb = 8.3,  SampleSizePre =  77, SampleSizePost =  69),
                               c(temp=23, muOb =  65.2, sigOb = 4.6,  SampleSizePre = 107, SampleSizePost =  96),
                               c(temp=27, muOb =  49.1, sigOb = 4.4,  SampleSizePre = 141, SampleSizePost = 114),
                               c(temp=30, muOb =  41.2, sigOb = 4.5,  SampleSizePre = 145, SampleSizePost = 131))

## Convert hours to days
MullensPupaeData           <- MullensPupaeDataHours
MullensPupaeData[,"muOb"]  <- MullensPupaeDataHours[,"muOb"]  / 24
MullensPupaeData[,"sigOb"] <- MullensPupaeDataHours[,"sigOb"] / 24
MullensPupaeData           <- as.data.frame(MullensPupaeData)

rm(MullensPupaeDataHours)

print("Mullens pupae data")
print(MullensPupaeData)
