# Flux Ratio - Active Index Surface Exchange (FRAISE)
# Thomas Loridan and Sue Grimmond
# ref: Loridan and Grimmond, 2011, JAMC.
#
# last modifications: TL 16/05/2011
#########################################################################################

#start with a clean work space
rm(list=ls())

#########################################################################################
# User input section 
# Default values are for the site of Lodz, Poland
#########################################################################################

# Path to the main directory (FRAISE) where the current script is located:
FRAISE_PATH <- "~/Energy/AustinAgua/UHI/FRAISE/"

# Output path (directory "FRAISE/outputs" by default):
PATH_OUT <- paste(FRAISE_PATH,"outputs/",sep="")

# Code used to name output files
SITE_CODE <- "LODZ"

# Day of year
DOY<- 172 # 21st of june

#General param 
LATITUDE <- 51.75
LONGITUDE <- 19.46
timezone <- 2 # including DST

# Urban param
URBAN_FRAC <- 0.7
ROOF_HEIGHT <- 10.6
ROOF_WIDTH <- 10.5
ROAD_WIDTH <- 14.1

#Vegetation param
GRASS_FRAC<-0.3
Spring_Start<-69
Spring_Stop<-144	
Fall_Start<-281
Fall_Stop<-324
LAI_GRASS <- 1.6

# LAI for trees
# Deciduous: (4->0) and coniferous (6->4) 
# here a mixture is assumed

LAI_MIN <-1
LAI_MAX <-6

# Source Active Surface Functions
source(paste(FRAISE_PATH,"util/Functions.R",sep=""))

NOON<-solar_noon(LATITUDE,LONGITUDE,timezone,DOY)
midday <- round(((NOON-floor(NOON))*24),0)

# Anthropogenic heat
# A mean midday estimate is needed here to predict flux values (QF_MEAN)
# In the example below it is taken from the LUCY model outputs (Allen et al, 2010, IJC)
# To simply provide an estimate comment the next four lines and assign a value to QF_MEAN
# If only interested in flux ratio values (and UZE) simply assign QF_MEAN<-0

QF_FILE <- paste(FRAISE_PATH,"data/QF_LUCY_",SITE_CODE,"_",DOY,".txt",sep="")
TABLE<-read.table(QF_FILE,skip=1,comment.char="%")
HOUR_LUCY<-TABLE[,2]
midday_indices <- which((HOUR_LUCY<=(midday+3))&(HOUR_LUCY>=(midday-3)))

QF_MEAN<-mean(TABLE[midday_indices,21])

# Estimate of mean midday incoming radiant energy (QDOWN_MEAN)
# Here it is estimated from observation data
# an alternative is to model both incoming shortwave and longwave radiation 
# from common meteorological data available from the US NCDC global database:
# http://cdo.ncdc.noaa.gov/pls/plclimprod/poemain.accessrouter?datasetabbv=DS3505&countryabbv=&georegionabbv=
# for simple models see for instance the following references:
# Kdown: Liston and Elder, 2006, Journal of Hydrometeorology
# Ldown: Loridan et al., 2011, JAMC. Offerle et al., 2003, JAMC
# 
# To simply provide an estimate comment the next six lines and assign a value to QDOWN_MEAN
# If only interested in flux ratio values (and UZE) simply assign QDOWN_MEAN<-0

QDOWN_FILE <- paste(FRAISE_PATH,"data/Qdown_",SITE_CODE,"_",DOY,".txt",sep="")
TABLE<-read.table(QDOWN_FILE,skip=0,comment.char="%")
HOUR_QDOWN<-TABLE[,6]
KDOWN <- TABLE[,15]
LDOWN <- TABLE[,18]
midday_indices <- which((HOUR_QDOWN<=(midday+3))&(HOUR_QDOWN>=(midday-3)))

QDOWN_MEAN<-mean(KDOWN[midday_indices]+LDOWN[midday_indices])

#########################################################################################
# End of user input section 
#########################################################################################

# FRAISE calculations
source(paste(FRAISE_PATH,"util/FRAISE_Calculations.R",sep=""))

# format and write output
OUT_DATA <- array(NA,dim=c(17,1))
OUT_DATA[1]  <- formatC(SITE_CODE,format="s",width=5)
OUT_DATA[2]  <- formatC(DOY,digits=3,format="d",width=5)
OUT_DATA[3]  <- formatC(UZE,digits=1,format="d",width=2)
OUT_DATA[4]  <- formatC(round(CHI_TOT,3),digits=3,format="f",width=6)
OUT_DATA[5]  <- formatC(round(CHI_URB,3),digits=3,format="f",width=6)
OUT_DATA[6]  <- formatC(round(CHI_VEG,3),digits=3,format="f",width=6)
OUT_DATA[7]  <- formatC(round(Qup_ratio,3),digits=3,format="f",width=6)
OUT_DATA[8]  <- formatC(round(Qs_ratio,3),digits=3,format="f",width=6)
OUT_DATA[9]  <- formatC(round(Qe_ratio,3),digits=3,format="f",width=6)
OUT_DATA[10] <- formatC(round(Bowen_ratio,3),digits=3,format="f",width=6)
OUT_DATA[11] <- formatC(round(QDOWN_MEAN,3),digits=3,format="f",width=8)
OUT_DATA[12] <- formatC(round(QF_MEAN,3),digits=3,format="f",width=8)
OUT_DATA[13] <- formatC(round(QUP_FRAISE,3),digits=3,format="f",width=8)
OUT_DATA[14] <- formatC(round(QS_FRAISE,3),digits=3,format="f",width=8)
OUT_DATA[15] <- formatC(round(QH_RESIDUAL,3),digits=3,format="f",width=8)
OUT_DATA[16] <- formatC(round(QH_BOWEN,3),digits=3,format="f",width=8)
OUT_DATA[17] <- formatC(round(QE_FRAISE,3),digits=3,format="f",width=8)
OUT_DATA[18] <- formatC(round(QSTAR_FRAISE,3),digits=3,format="f",width=8)
			
TXT_PATH <- paste(PATH_OUT,SITE_CODE,"_",DOY,"_FRAISE_OUTPUT.txt",sep="")
if ((file.exists(TXT_PATH)==FALSE)){
	Output_Header(TXT_PATH)
	write(OUT_DATA,file=TXT_PATH,ncolumns=length(OUT_DATA),append = TRUE, sep = "   ")
}else {
	print(paste("File: ",TXT_PATH," already exists (delete first if you want to overwrite)",sep=""))
}

# End of program



