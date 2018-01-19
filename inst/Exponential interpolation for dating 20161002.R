### PLEASE READ:

# This is a routine used to fill-in gaps for Pb-210 dating. 
# Columns must contain, in this order: Sample code, layer depth = section bottom depth (cm), 
# its uncertainty (e.g. ruler resolution), section mass (g), its uncertainty (balance resolution),
# Pb-210 (Bq kg-1), its uncertainty, Ra-226, its uncertainty.


# If blanks are found at the start or end, the code will not fill them. 
# Values at the end might not be needed, but values at the start are compulsory for CF dating.
# The user should decide how to solve that serious problem. 

# For the meaning of variables, and to cite this work (please!):
# Sanchez-Cabeza, J. A., Ruiz-Fernández, A. C., Ontiveros-Cuadras, J. F., Bernal, L. H. P., Olid, C. 
# (2014). Monte Carlo uncertainty calculation of 210 Pb chronologies and accumulation rates of 
# sediments and peat bogs. Quaternary Geochronology, 23, 80-93.

# Note: nothing is better than a complete analysis of a core 
# and a careful check of all results before attempting any dating. 

### START

rm(list=ls()) # Clear all objects
library(zoo) # Used to fill in random gaps in all variables (except mass), cumulative sums and means
library(beepr)

### DATA TO BE REVISED BY THE USER
setwd("C:/Users/unam/Documents/RStudio/Dating Pb210/MINECO/ABRA3")
# Remember to change backslah by /
NameCore="ABRA3 20171006"

### File names
NameInput=paste0(NameCore,".csv") #paste0 does not add spaces
NameOutput=paste(NameCore,"complete.csv")

### Reads data
Table=read.csv(NameInput,header=F,skip=1,stringsAsFactors=F,check.names=F) # Table with header
Variables=c("SampleCode","DepthLayer","DepthLayerUnc","MassSection","MassSectionUnc",
	"Pb210","Pb210Unc","Ra226","Ra226Unc") # To make a “nice” Table
colnames(Table)=Variables
SampleCode=Table$SampleCode
MassSection=Table$MassSection

# Use the section accumulated mass for interpolation (core surface not needed)
MassLayerAcc=c(0,cumsum(MassSection))
MassSectionAcc=rollmean(MassLayerAcc,2) # Rolling mean

### Create a zoo series and make exponential interpolation
TableZoo=read.zoo(data.frame(cbind(MassSectionAcc,Table[,2:9])),index.column = 1)
TableZooLog=log(TableZoo) # For exponential interpolation
TableZooLogFull=na.approx(TableZooLog) # Fill gaps
TableZooFull=exp(TableZooLogFull) # Return to linear

### Save results
write.table(cbind(SampleCode,coredata(TableZooFull)),
       file = NameOutput,sep=",",append=F,row.names=F,quote=F)
beep(1) 
### END


