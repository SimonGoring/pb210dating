#### PLEASE READ: ####
# This is a routine designed to provide CF and CFCS Pb-210 chronologies.
# Columns must contain, in this order:
#    Sample code, layer depth = section bottom depth (cm), its uncertainty (e.g. ruler resolution),
#    section mass (g), its uncertainty (balance resolution), Pb-210 (Bq kg-1), its uncertainty,
#    Ra-226 (Bq kg-1), its uncertainty.
# All sections used for dating must have complete information (no gaps allowed), including Ra-226.
# For the meaning of variables, and to cite this work (please!):
# 1. Sanchez-Cabeza, J. A., & Ruiz-Fern?ndez, A. C. (2012).
# 210 Pb sediment radiochronology: an integrated formulation and classification of dating models.
# Geochimica et Cosmochimica Acta, 82, 183-200.
# 2. Sanchez-Cabeza, J. A., Ruiz-Fern?ndez, A. C., Ontiveros-Cuadras, J. F., Bernal, L. H. P., Olid, C. (2014).
# Monte Carlo uncertainty calculation of 210 Pb chronologies and accumulation rates of # sediments and peat bogs.
# Quaternary Geochronology, 23, 80-93.
#
# Note: nothing is better than the complete analysis of a core
# and a careful check of all results before attempting any dating.

#### START ####
rm(list=ls()) # Clear all objects
library(zoo) # Used for cumulative sums and means
library(beepr)
library(lubridate)

#### DATA TO BE REVISED BY THE USER ####

Scientist="JASC" # Enter your name or initials
setwd("C:/Users/unam/Documents/RStudio/Dating Pb210/MINECO/ABRA2") # Enter your working directory
# Remember to change backslash for /
NameCore="ABRA2 20171006 v2 complete" # Keep the quotes
TimeSampling=as.Date("2015-09-15") # Enter yyyy-mm-dd, keep the quotes.
Diameter=8.5 # cm
DiameterUnc=0.1/2/sqrt(3) # Instrumental resolution
NSectionsDating=38 # Number of sections to be used for dating.
MissingInventoryTop=0.0 # In Bq. Usually zero, or estimated using a spreadsheet.
MissingInventoryBottom=7.46 # In Bq. Usually zero, or estimated using a spreadsheet.
Runs=100000 # At this moment, 100.000 is recommended.

#### Constants ####

TimeStart=proc.time()
Halflife=22.23; HalflifeUnc=0.12 # DDEP (2016)
ReductionFactor=100 # Can try others, this one works
YearSampling=decimal_date(TimeSampling)
AgeSampling=as.Date(TimeSampling) # For final output

#### File names ####

NameInput=paste0(NameCore,".csv") # paste0 does not add spaces
NameDatingCF=paste(NameCore,"Dating CF.csv")
NameDatingCF2=paste(NameCore,"Dating CF user.csv")
NameDatingCFCS=paste(NameCore,"Dating CFCS.csv")
NameMetadata=paste(NameCore,"Metadata.csv")
# NamePlot=paste(NameCore,"Plots.pdf") # Pending

#### Reads data and vectors ####

Table=read.csv(NameInput,header=F,skip=1,stringsAsFactors=F,check.names=F) # Table with header, not used
NSections=length(Table$V1) # Total number of sections
SampleCode=Table$V1
DepthLayer=c(0,Table$V2)
DepthLayerUnc=c(0,Table$V3)
DepthLayerUnc=DepthLayerUnc/sqrt(3) # Square distribution (optional line)
MassSection=Table$V4
MassSectionUnc=Table$V5
MassSectionUnc=MassSectionUnc/sqrt(3) # Square distribution (optional line)
Pb210=Table$V6
Pb210Unc=Table$V7
Ra226=Table$V8
Ra226Unc=Table$V9
Variables=c("SampleCode","DepthLayer","DepthLayerUnc","MassSection","MassSectionUnc",
            "Pb210","Pb210Unc","Ra226","Ra226Unc") # For the internal Table use
colnames(Table)=Variables

#### Dating calculations ####

Lambda=log(2)/Halflife
DiameterUnc=DiameterUnc/sqrt(3) # Square distribution (optional line)
Surface=pi*(Diameter/2)^2 # Change if the core is not a cylinder
DepthDelta=diff(DepthLayer)
DepthSection=rollmean(DepthLayer,2) # Rolling mean
Density=MassSection/Surface/DepthDelta
MassDepthSection=MassSection/Surface
MassDepthLayerAcc=c(0,cumsum(MassDepthSection)) # c(0,...) is to complete the surface layer
MassDepthSectionAcc=rollmean(MassDepthLayerAcc,2)
Concentration=Pb210[1:NSectionsDating]-Ra226[1:NSectionsDating] # Only for dating sections

Deposit=Concentration*MassSection[1:NSectionsDating]/1000
Deposit=c(Deposit,MissingInventoryBottom) # Adds the bottom missing inventory, as an extra section
DepositAbove=cumsum(c(MissingInventoryTop,Deposit)) # Adds the top missing inventory to the first layer
DepositBelow=c(rev(cumsum(rev(Deposit))),0) # Calculates accumulated inventory in reverse order
DepositBelow[1]=DepositBelow[1]+MissingInventoryTop # Adds top missing inventory to first layer
#DepositAbove[length(DepositAbove)]==DepositBelow[1] # Test inventories are equal (TRUE)
Inventory=DepositBelow[1]

Flux=Lambda*Inventory/Surface*10000
TimeLayer=log(1+DepositAbove/DepositBelow)/Lambda
length(TimeLayer)=NSectionsDating+1 # Last layer has no date (infinity)
TimeSection=rollmean(TimeLayer,2)
TimeDelta=diff(TimeLayer)
MarSection=MassSection[1:NSectionsDating]/TimeDelta/Surface # Last section cannot be dated
SarSection=DepthDelta[1:NSectionsDating]/TimeDelta

#### CFCS calculation of mean MAR and SAR ####

ConcentrationCF=Concentration
ConcentrationCF[Concentration==0]=NA
MARRegression=lm(log(ConcentrationCF)~MassDepthSectionAcc[1:NSectionsDating],na.action=na.exclude)
MARRegressionSummary=summary(MARRegression)
MARSlope=MARRegressionSummary$coefficients[2,1]
MARSlopeUnc=MARRegressionSummary$coefficients[2,2]
MARSlopeT=MARRegressionSummary$coefficients[2,3]
MARSlopeP=MARRegressionSummary$coefficients[2,4]

SARRegression=lm(log(ConcentrationCF)~DepthSection[1:NSectionsDating],na.action=na.exclude)
SARRegressionSummary=summary(SARRegression)
SARSlope=SARRegressionSummary$coefficients[2,1]
SARSlopeUnc=SARRegressionSummary$coefficients[2,2]
SARSlopeT=SARRegressionSummary$coefficients[2,3]
SARSlopeP=SARRegressionSummary$coefficients[2,4]

MARMean=-Lambda/MARSlope
SARMean=-Lambda/SARSlope
TimeCFCS=MassDepthSectionAcc/MARMean

#### Clear variances ####

SurfaceVar=0
LambdaVar=0
DepthDeltaVar=0
DepthSectionVar=0
DensityVar=0
MassDepthSectionVar=0
MassDepthLayerAccVar=0
MassDepthSectionAccVar=0
ConcentrationVar=0
DepositVar=0
DepositAboveVar=0
DepositBelowVar=0
FluxVar=0
TimeLayerVar=0
TimeSectionVar=0
TimeDeltaVar=0
MarSectionVar=0
SarSectionVar=0
MARMeanVar=0
SARMeanVar=0
TimeCFCSVar=0

#### Monte Carlo calculations: ####

for(i in 1:Runs) {

  # Generation of random values
  	# For measured uncertainty, generate random number from normal distribution
  DiameterRan=rnorm(1,Diameter,DiameterUnc/ReductionFactor) # Reduction factor avoids divergences
  	# For calculated uncertainty, calculation is repeated with the random numbers
  SurfaceRan=pi*(DiameterRan/2)^2
  HalflifeRan=rnorm(1,Halflife,HalflifeUnc/ReductionFactor)
  LambdaRan=log(2)/HalflifeRan
  DepthLayerRan=rnorm(NSections+1,DepthLayer,DepthLayerUnc/ReductionFactor) # 1 more layer than sections
  DepthDeltaRan=diff(DepthLayerRan)
  DepthSectionRan=rollmean(DepthLayerRan,2)
  MassSectionRan=rnorm(NSections,MassSection,MassSectionUnc/ReductionFactor)
  DensityRan=MassSectionRan/SurfaceRan/DepthDeltaRan
  MassDepthSectionRan=MassSectionRan/SurfaceRan
  MassDepthLayerAccRan=c(0,cumsum(MassDepthSectionRan))
  MassDepthSectionAccRan=rollmean(MassDepthLayerAccRan,2)
  Pb210Ran=rnorm(NSections,Pb210,Pb210Unc/ReductionFactor)
  Ra226Ran=rnorm(NSections,Ra226,Ra226Unc/ReductionFactor)
  ConcentrationRan=Pb210Ran[1:NSectionsDating]-Ra226Ran[1:NSectionsDating]
  DepositRan=ConcentrationRan*MassSectionRan[1:NSectionsDating]/1000
  DepositRan=c(DepositRan,MissingInventoryBottom)
  DepositAboveRan=cumsum(c(MissingInventoryTop,DepositRan))
  DepositBelowRan=c(rev(cumsum(rev(DepositRan))),0)
  DepositBelowRan[1]=DepositBelowRan[1]+MissingInventoryTop
  InventoryRan=DepositBelowRan[1]
  FluxRan=LambdaRan*InventoryRan/SurfaceRan*10000
  TimeLayerRan=log(1+DepositAboveRan/DepositBelowRan)/LambdaRan
  length(TimeLayerRan)=NSectionsDating+1 # Last layer has no date (infinity)
  TimeSectionRan=rollmean(TimeLayerRan,2)
  TimeDeltaRan=diff(TimeLayerRan)
  MarSectionRan=MassSectionRan[1:NSectionsDating]/TimeDeltaRan/SurfaceRan
  SarSectionRan=DepthDeltaRan[1:NSectionsDating]/TimeDeltaRan
  MARMeanRan=-LambdaRan/rnorm(1,MARSlope,MARSlopeUnc/ReductionFactor)
  SARMeanRan=-LambdaRan/rnorm(1,SARSlope,SARSlopeUnc/ReductionFactor)
  TimeCFCSRan=MassDepthSectionAccRan/MARMeanRan

#### Caculation of variances ####

  SurfaceVar=SurfaceVar+(Surface-SurfaceRan)^2
  LambdaVar=LambdaVar+(Lambda-LambdaRan)^2
  DepthDeltaVar=DepthDeltaVar+(DepthDelta-DepthDeltaRan)^2
  DepthSectionVar=DepthSectionVar+(DepthSection-DepthSectionRan)^2
  DensityVar=DensityVar+(Density-DensityRan)^2
  MassDepthSectionVar=MassDepthSectionVar+(MassDepthSection-MassDepthSectionRan)^2
  MassDepthLayerAccVar=MassDepthLayerAccVar+(MassDepthLayerAcc-MassDepthLayerAccRan)^2
  MassDepthSectionAccVar=MassDepthSectionAccVar+(MassDepthSectionAcc-MassDepthSectionAccRan)^2
  ConcentrationVar=ConcentrationVar+(Concentration-ConcentrationRan)^2
  DepositVar=DepositVar+(Deposit-DepositRan)^2
  DepositAboveVar=DepositAboveVar+(DepositAbove-DepositAboveRan)^2
  DepositBelowVar=DepositBelowVar+(DepositBelow-DepositBelowRan)^2
  FluxVar=FluxVar+(Flux-FluxRan)^2
  TimeLayerVar=TimeLayerVar+(TimeLayer-TimeLayerRan)^2
  # length(TimeLayerVar)=NSectionsDating+1 # Last layer has no date (infinity). Not needed?
  TimeSectionVar=TimeSectionVar+(TimeSection-TimeSectionRan)^2
  TimeDeltaVar=TimeDeltaVar+(TimeDelta-TimeDeltaRan)^2
  MarSectionVar=MarSectionVar+(MarSection-MarSectionRan)^2
  SarSectionVar=SarSectionVar+(SarSection-SarSectionRan)^2
  MARMeanVar =MARMeanVar+(MARMean-MARMeanRan)^2
  SARMeanVar =SARMeanVar+(SARMean-SARMeanRan)^2
  TimeCFCSVar=TimeCFCSVar+(TimeCFCS-TimeCFCSRan)^2
  }

#### Calculate final uncertainties ####

SurfaceUnc=sqrt(SurfaceVar/Runs)*ReductionFactor
LambdaUnc=sqrt(LambdaVar/Runs)*ReductionFactor
DepthDeltaUnc=sqrt(DepthDeltaVar/Runs)*ReductionFactor
DepthSectionUnc=sqrt(DepthSectionVar/Runs)*ReductionFactor
DensityUnc=sqrt(DensityVar/Runs)*ReductionFactor
MassDepthSectionUnc=sqrt(MassDepthSectionVar/Runs)*ReductionFactor
MassDepthLayerAccUnc=sqrt(MassDepthLayerAccVar/Runs)*ReductionFactor
MassDepthSectionAccUnc=sqrt(MassDepthSectionAccVar/Runs)*ReductionFactor
ConcentrationUnc=sqrt(ConcentrationVar/Runs)*ReductionFactor
DepositUnc=sqrt(DepositVar/Runs)*ReductionFactor
DepositAboveUnc=sqrt(DepositAboveVar/Runs)*ReductionFactor
DepositBelowUnc=sqrt(DepositBelowVar/Runs)*ReductionFactor
FluxUnc=sqrt(FluxVar/Runs)*ReductionFactor
TimeLayerUnc=sqrt(TimeLayerVar/Runs)*ReductionFactor
# length(TimeLayerUnc)=NSectionsDating+1 # Last layer has no date (infinity). Not needed?
TimeSectionUnc=sqrt(TimeSectionVar/Runs)*ReductionFactor
TimeDeltaUnc=sqrt(TimeDeltaVar/Runs)*ReductionFactor
MarSectionUnc=sqrt(MarSectionVar/Runs)*ReductionFactor
SarSectionUnc=sqrt(SarSectionVar/Runs)*ReductionFactor
InventoryUnc=DepositBelowUnc[1]
MARMeanUnc=sqrt(MARMeanVar/Runs)*ReductionFactor
SARMeanUnc=sqrt(SARMeanVar/Runs)*ReductionFactor
TimeCFCSUnc=sqrt(TimeCFCSVar/Runs)*ReductionFactor

#### Presentation of results ####

#### Metadata file ####

# Elapsed time
TimeEnd=proc.time()
TimeElapsed=TimeEnd[3]-TimeStart[3]

MetadataRows=c("Scientist","Date","Core name","Sampling Date","Diameter(cm)","u(Diameter)",
               "Number of Sections","Dated Sections (CF)","Half life (yr)","uHalf life (yr)",
               "Missing inventory top (%)","Missing inventory bottom (%)",
               "210Pb Flux (Bq m-2 yr-1)","u(210Pb flux)",
               "Mean MAR (g cm-2 yr-1)","u(mean MAR)","Mean SAR (cm yr-1)","u(mean SAR)",
               "Runs","Time elapsed (s)")
Metadata=c(Scientist,as.character(Sys.time()),NameCore,as.character(AgeSampling),Diameter,DiameterUnc,
           as.integer(NSections),as.integer(NSectionsDating),Halflife,HalflifeUnc,
           MissingInventoryTop/Inventory*100,MissingInventoryBottom/Inventory*100,
           Flux,FluxUnc,
           MARMean,MARMeanUnc,SARMean,SARMeanUnc,
           as.integer(Runs),TimeElapsed)
write.table(Metadata,sep=",",file=NameMetadata,append=F,row.names=MetadataRows,col.names=F,quote=F)

#### CF Dating complete file ####

YearLayer=YearSampling-TimeLayer
AgeLayer=as.Date(TimeSampling-TimeLayer*365.24) # Calculated
TimeTop=TimeLayer[1:NSectionsDating]
YearTop=YearLayer[1:NSectionsDating]
AgeTop=AgeLayer[1:NSectionsDating]
AgeTopUnc=TimeLayerUnc[1:NSectionsDating]
TimeBottom=TimeLayer[-1]
YearBottom=YearLayer[-1]
AgeBottom=AgeLayer[-1]
AgeBottomUnc=TimeLayerUnc[-1]
YearSection=YearSampling-TimeSection
AgeSection=as.Date(TimeSampling-TimeSection*365.25) # Calculated
AgeSectionUnc=TimeSectionUnc

# Data.Frames are rectangular. Adjust length to NSections. Note: length transforms Dates to numbers
NSectionsTable=ifelse(MissingInventoryBottom==0,NSectionsDating-1,NSectionsDating)
CompleteTable=function(Column) {              # Function to clean the table (completes with NA)
  Column[(NSectionsTable+1):NSections]=NA       # Deletes from the next section
  Column
}
if (NSections!=NSectionsDating) {
  Concentration[(NSectionsDating+1):NSections]=NA
  ConcentrationUnc[(NSectionsDating+1):NSections]=NA
  TimeTop=CompleteTable(TimeTop)
  YearTop=CompleteTable(YearTop)
  AgeTop=CompleteTable(AgeTop)
  AgeTopUnc=CompleteTable(AgeTopUnc)
  TimeSection=CompleteTable(TimeSection)
  YearSection=CompleteTable(YearSection)
  AgeSection=CompleteTable(AgeSection)
  AgeSectionUnc=CompleteTable(AgeSectionUnc)
  TimeBottom=CompleteTable(TimeBottom)
  YearBottom=CompleteTable(YearBottom)
  AgeBottom=CompleteTable(AgeBottom)
  AgeBottomUnc=CompleteTable(AgeBottomUnc)
  MarSection=CompleteTable(MarSection)
  MarSectionUnc=CompleteTable(MarSectionUnc)
  SarSection=CompleteTable(SarSection)
  SarSectionUnc=CompleteTable(SarSectionUnc)
}

# Table heading
DatingColumnsCF=c("Sample code","Depth (cm)","Mass depth (g/cm2)","Density (g/cm3)","+-",
                "210Pb (Bq yr-1)","+-","226Ra (Bq yr-1)","+-","210Pbex (Bq yr-1)","+-",
                "Top time (yr)","Top time (yr)","Top age (A.D)","+-",
                "Mean time (yr)","Mean time (yr)","Mean age (A.D.)","+-",
                "Bottom time (yr)","Bottom time (yr)","Bottom age (A.D)","+-",
                "MAR (g/cm2 yr)","+-","SAR (cm/yr)","+-")

TableDatingCF=data.frame(SampleCode,DepthSection,MassDepthSectionAcc,Density,DensityUnc,
  Pb210,Pb210Unc,Ra226,Ra226Unc,Concentration,ConcentrationUnc,
  TimeTop,YearTop,AgeTop,AgeTopUnc,TimeSection,YearSection,AgeSection,AgeSectionUnc,
  TimeBottom,YearBottom,AgeBottom,AgeBottomUnc,
  MarSection,MarSectionUnc,SarSection,SarSectionUnc)

write.table(as.data.frame(TableDatingCF),file=NameDatingCF,sep=",",append=F,row.names=F,
            col.names=DatingColumnsCF,quote=F,na="")

#### CF Dating user file ####

DatingColumnsCF2=c("Sample code","Depth (cm)",
                  "Top time (yr)","+-",
                  "Mean time (yr)","+-",
                  "Bottom time (yr)","+-",
                  "MAR (g/cm2 yr)","+-",
                  "SAR (cm/yr)","+-")

TableDatingCF2=data.frame(SampleCode,DepthSection,YearTop,AgeTopUnc,YearSection,AgeSectionUnc,
                         YearBottom,AgeBottomUnc,MarSection,MarSectionUnc,SarSection,SarSectionUnc)

write.table(as.data.frame(TableDatingCF2),file=NameDatingCF2,sep=",",append=F,row.names=F,
            col.names=DatingColumnsCF2,quote=F,na="")

#### CFCS Dating file ####

YearCFCS=YearSampling-TimeCFCS
AgeCFCS=as.character(TimeSampling-TimeCFCS*365.25) # Calculated
AgeCFCSUnc=TimeCFCSUnc

# Table heading
DatingColumnsCFCS=c("Sample code","Depth (cm)","Mass depth (g/cm2)","Density (g/cm3)","+-",
                  "210Pb (Bq yr-1)","+-","226Ra (Bq yr-1)","+-","210Pbex (Bq yr-1)","+-",
                  "Time CFCS (yr)","Year CFCS (yr)","Age CFCS (A.D)","+-","MAR (g/cm2 yr","+-","SAR (cm/yr)","+-")

# Equal length needed for write.table
length(SampleCode)=length(DepthSection)=length(MassDepthSectionAcc)=length(Density)=length(DensityUnc)=
  length(Pb210)=length(Pb210Unc)=length(Ra226)=length(Ra226Unc)=length(Concentration)=length(ConcentrationUnc)=
  length(TimeCFCS)=length(YearCFCS)=length(AgeCFCS)=length(AgeCFCSUnc)=length(MARMean)=length(MARMeanUnc)=
  length(SARMean)=length(SARMeanUnc)=NSections

TableDatingCFCS=cbind(SampleCode,DepthSection,MassDepthSectionAcc,Density,DensityUnc,
                    Pb210,Pb210Unc,Ra226,Ra226Unc,Concentration,ConcentrationUnc,
                    TimeCFCS,YearCFCS,AgeCFCS,AgeCFCSUnc,MARMean,MARMeanUnc,SARMean,SARMeanUnc)
# NA are set to blanks
TableDatingCFCS[is.na(TableDatingCFCS)]=""
# TableDating[is.nan(TableDating)]="" Does not work
write.table(TableDatingCFCS,file=NameDatingCFCS,sep=",",append=F,row.names=F,col.names=DatingColumnsCFCS,quote=F)

#### END ####
beep(3)
