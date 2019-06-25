# Program called by main.R

# Normalized roof height
h<-ROOF_HEIGHT/(ROOF_WIDTH+ROAD_WIDTH)

# Normalized roof width
r<-ROOF_WIDTH/(ROOF_WIDTH+ROAD_WIDTH)

# Normalized road width
w<-ROAD_WIDTH/(ROOF_WIDTH+ROAD_WIDTH)

lst_time<-seq(1,24)

# Active portion of surfaces (shading patterns etc...)
ActiveSurfaces <- Active_Surface_Urban(h,r,w,LATITUDE,LONGITUDE,timezone,DOY)
active_road<-as.numeric(ActiveSurfaces[1,])
active_wall<-as.numeric(ActiveSurfaces[2,])
active_roof<-as.numeric(ActiveSurfaces[3,])

# Vegetation seasonal evolution
V_phenology<-Vegetation_Phenology(Spring_Start,Spring_Stop,Fall_Start,Fall_Stop,LATITUDE,DOY)
LAI_TREE<-(V_phenology*(LAI_MAX-LAI_MIN))+LAI_MIN

# Link to the overall 3d surface
ACTIVE_INDEX<-Active_Surface_Index(active_road,active_wall,active_roof,
URBAN_FRAC,ROOF_HEIGHT,ROAD_WIDTH,ROOF_WIDTH,GRASS_FRAC,LAI_TREE,LAI_MAX,LAI_GRASS)

# Get active surface indices
CHI_URB<-mean(as.numeric(ACTIVE_INDEX[1,(midday-3):(midday+3)]))
CHI_VEG<-mean(as.numeric(ACTIVE_INDEX[2,(midday-3):(midday+3)]))
CHI_TOT<-mean(as.numeric(ACTIVE_INDEX[3,(midday-3):(midday+3)]))

# Empirical coefficient for predictive Flux Ratio relations
A   <- list(Qup= 1.04, Qs=-1.90, Qe=-0.200, Bowen= 3.80)
B   <- list(Qup=-0.15, Qs= 0.89, Qe= 0.016, Bowen=-0.22)
C   <- list(Qup= 0.62, Qs= 0.13, Qe= 0.110, Bowen= 1.60)
X0  <- list(Qup= 0.35, Qs= 0.11, Qe= 0.430, Bowen= 0.43)
#
CHI <- list(Qup= CHI_TOT, Qs= CHI_URB, Qe= CHI_VEG, Bowen= CHI_VEG)

Qup_ratio <- A$Qup*pmax(0,(X0$Qup-CHI$Qup))+B$Qup*pmax(0,(CHI$Qup-X0$Qup))+C$Qup
Qs_ratio <- A$Qs*pmax(0,(X0$Qs-CHI$Qs))+B$Qs*pmax(0,(CHI$Qs-X0$Qs))+C$Qs
Qe_ratio <- A$Qe*pmax(0,(X0$Qe-CHI$Qe))+B$Qe*pmax(0,(CHI$Qe-X0$Qe))+C$Qe
Bowen_ratio <-A$Bowen*pmax(0,(X0$Bowen-CHI$Bowen))+B$Bowen*pmax(0,(CHI$Bowen-X0$Bowen))+C$Bowen
#
Qh_ratio <- 1 - Qup_ratio - Qs_ratio - Qe_ratio

# Compute fluxes and add anthropogenic heat contribution
QF_COEFF <- list(alpha=0.2, beta=0.6, gamma=0.2)

QUP_FRAISE <- Qup_ratio*QDOWN_MEAN+QF_COEFF$beta*QF_MEAN
QS_FRAISE <- Qs_ratio*QDOWN_MEAN+QF_COEFF$gamma*QF_MEAN
QE_FRAISE <- Qe_ratio*QDOWN_MEAN
QSTAR_FRAISE <- QDOWN_MEAN - QUP_FRAISE

QH_RESIDUAL <- Qh_ratio*QDOWN_MEAN +QF_COEFF$alpha*QF_MEAN 
QH_BOWEN <- (QE_FRAISE*Bowen_ratio)+QF_COEFF$alpha*QF_MEAN

# UZE

CHI_THRESH <- list(URB= 0.11, VEG=0.43)

if((CHI_URB <  CHI_THRESH$URB)&(CHI_VEG <  CHI_THRESH$VEG)) { UZE <- 1 } 
if((CHI_URB <  CHI_THRESH$URB)&(CHI_VEG >= CHI_THRESH$VEG)) { UZE <- 2 } 
if((CHI_URB >= CHI_THRESH$URB)&(CHI_VEG >= CHI_THRESH$VEG)) { UZE <- 3 } 
if((CHI_URB >= CHI_THRESH$URB)&(CHI_VEG <  CHI_THRESH$VEG)) { UZE <- 4 } 




