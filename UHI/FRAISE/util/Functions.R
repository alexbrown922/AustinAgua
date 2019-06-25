
solar_noon<- function(LATITUDE,LONGITUDE,timezone,DOY) {
			lst_time<-seq(1,1440)
			DECTIME<-DOY+lst_time/1440.

			# Oaf - NARP
			DEG2RAD<-0.017453292
			RAD2DEG<-57.29577951
			latr<-DEG2RAD*LATITUDE
			lonr<-DEG2RAD*LONGITUDE
			
			#compute zenith angle (from NARP)
			gamma<-2*3.141592654/365.25463*DECTIME
			eqtime<-229.18*(7.5e-5+1.868e-3*cos(gamma)-0.032077*sin(gamma)-0.014615*cos(2.*gamma)-0.040849*sin(2.*gamma))
			delta<-6.918e-3-0.399912*cos(gamma)+0.070257*sin(gamma)-0.006758*cos(2.*gamma)+9.07e-4*sin(2.*gamma)-2.697e-3*cos(3.*gamma)+1.48e-3*sin(3.*gamma)

			#Oaf-NARP
			# Oke (p340): for local mean solar time, ADD 4 min for each degree EAST of longitude
			#Then add the eqtime to get local apparent solar time
			time_offset<-eqtime+4.*lonr*RAD2DEG-60.*timezone

			tst<-(DECTIME-floor(DECTIME))*1440.+time_offset
			ha<- 15*(12-tst/60)
			wt<-ha*DEG2RAD	   

			THETA_Z<-abs(acos(sin(latr)*sin(delta)+cos(latr)*cos(delta)*cos(wt)))
			NOON<-DECTIME[which(THETA_Z==min(THETA_Z))]
			return(NOON)
}

toolkit_FilterData <- function(TABLE,REF_COL,THRES_MIN,THRES_MAX) {
	# REF_COL, the column of the matrix on which we base the filtering criteria
	# removes the data outside of the window [THRES_MIN -> THRES_MAX]
	# INPUT: TABLE, the matrix to filter
	# THRES_MIN,THRES_MAX, the thresholds to use
	# TXT_FILE for outputing text information
	################################################
	REF_TMP <- TABLE[,REF_COL]
	NROW <- length(REF_TMP)
	NCOLUMN <- length(TABLE[1,])
	REF_TMP[(REF_TMP < THRES_MIN | REF_TMP > THRES_MAX)] <- NA
	TABLE_TMP <- TABLE[which(is.finite(REF_TMP)),]
	return(TABLE_TMP)
}


Active_Surface_Urban <- function(h,r,w,LATITUDE,LONGITUDE,timezone,DOY) {	
			lst_time<-seq(1,24)
			DECTIME<-DOY+lst_time/24

			# Oaf - NARP
			DEG2RAD<-0.017453292
			RAD2DEG<-57.29577951
			latr<-DEG2RAD*LATITUDE
			lonr<-DEG2RAD*LONGITUDE

			#compute zenith angle (from NARP)
			gamma<-2*3.141592654/365.25463*DECTIME
			eqtime<-229.18*(7.5e-5+1.868e-3*cos(gamma)-0.032077*sin(gamma)-0.014615*cos(2.*gamma)-0.040849*sin(2.*gamma))
			delta<-6.918e-3-0.399912*cos(gamma)+0.070257*sin(gamma)-0.006758*cos(2.*gamma)+9.07e-4*sin(2.*gamma)-2.697e-3*cos(3.*gamma)+1.48e-3*sin(3.*gamma)

			#Oaf-NARP
			# Oke (p340): for local mean solar time, ADD 4 min for each degree EAST of longitude
			#Then add the eqtime to get local apparent solar time
			time_offset<-eqtime+4.*lonr*RAD2DEG-60.*timezone

			tst<-(DECTIME-DOY)*1440.+time_offset
			ha<- 15*(12-tst/60)
			wt<-ha*DEG2RAD	   

			THETA_Z<-abs(acos(sin(latr)*sin(delta)+cos(latr)*cos(delta)*cos(wt)))
			THETA_S<-abs(asin(sin(latr)*sin(delta)+cos(latr)*cos(delta)*cos(wt)))

			# Also find the value of wt at sunrise and sunset (i.e THETA_Z=pi/2)
			wt_ss=acos((cos(pi/2)-sin(latr)*sin(delta))/(cos(latr)*cos(delta)))
			hss=wt_ss*RAD2DEG
			day_hours<-2*mean(hss)/15
			day_hours_index<-which(THETA_Z<pi/2)
			day_night_index<-which(THETA_Z>=pi/2)

			##############################################
			#Kusaka's SLUCM (see paper and code)

			PHI<-LATITUDE*DEG2RAD

			SNT<-(cos(delta)*sin(wt))/cos(THETA_S)
			CNT<-(cos(THETA_Z)*sin(PHI)-sin(delta))/cos(THETA_S)/cos(PHI)

			#Consider 8 different canyon orientations
			HOUI<-array(NA,dim=c(length(SNT),8))

			HOUI[,1]<-(SNT*cos(pi/8.)   -CNT*sin(pi/8.))
			HOUI[,2]<-(SNT*cos(2.*pi/8.)-CNT*sin(2.*pi/8.))
			HOUI[,3]<-(SNT*cos(3.*pi/8.)-CNT*sin(3.*pi/8.))
			HOUI[,4]<-(SNT*cos(4.*pi/8.)-CNT*sin(4.*pi/8.))
			HOUI[,5]<-(SNT*cos(5.*pi/8.)-CNT*sin(5.*pi/8.))
			HOUI[,6]<-(SNT*cos(6.*pi/8.)-CNT*sin(6.*pi/8.))
			HOUI[,7]<-(SNT*cos(7.*pi/8.)-CNT*sin(7.*pi/8.))
			HOUI[,8]<-(SNT*cos(8.*pi/8.)-CNT*sin(8.*pi/8.))
			
			# Active road 
			AR<-w-(h*abs(tan(THETA_Z))*abs(HOUI))
			AR[which(AR<0)]<-0
			AR[day_night_index,]<-0

			active_road<-(AR[,1]+AR[,2]+AR[,3]+AR[,4]+AR[,5]+AR[,6]+AR[,7]+AR[,8])/8.

			# Active wall
			THETA_ZLIM<-atan(w/h)
			THETA_Z_WALL<-THETA_Z
			THETA_Z_WALL[which(THETA_Z_WALL<=THETA_ZLIM)]<-THETA_ZLIM
			AW<-w/abs(tan(THETA_Z_WALL))*abs(HOUI)
			AW[which(AW<0)]<-0
			AW[day_night_index,]<-0

			active_wall<-(AW[,1]+AW[,2]+AW[,3]+AW[,4]+AW[,5]+AW[,6]+AW[,7]+AW[,8])/8.

			AROOF<-0*AW
			AROOF[which(AW>0)]<-r
			AROOF[which(AR>0)]<-r
			AROOF[day_night_index,]<-0
			AROOF[which(AROOF<0)]<-0

			active_roof<-(AROOF[,1]+AROOF[,2]+AROOF[,3]+AROOF[,4]+AROOF[,5]+AROOF[,6]+AROOF[,7]+AROOF[,8])/8.

			return(rbind(active_road,active_wall,active_roof))
}#end function

Vegetation_Phenology <- function(Spring_Start,Spring_Stop,Fall_Start,Fall_Stop,LATITUDE,DOY) {
		
	
		M0<-0.03
		ds<-(Spring_Start+Spring_Stop)/2
		df<-(Fall_Start+Fall_Stop)/2
		kks=log10((1-M0)/M0)/(ds-Spring_Start) #Spring
		kkf=log10((1-M0)/M0)/(Fall_Stop-df) #Fall
	
		if(LATITUDE>0) {
			V_phenology<-(1/(1+10**(kks*(ds-DOY))))*(1/(1+10**(kkf*(DOY-df))))
		}
	
		if(LATITUDE<0) {
			V_phenology<-((1/(1+10**(kks*(ds-DOY))))+(1/(1+10**(kkf*(DOY-df)))))
		}
	
		return(V_phenology)

}#end function

Active_Surface_Index <- function(active_road,active_wall,active_roof,URBAN_FRAC,ROOF_HEIGHT,ROAD_WIDTH,ROOF_WIDTH,GRASS_FRAC,LAI_TREE,LAI_MAX,LAI_GRASS) {

		# In terms of material surfaces sun lit ...
		# The complete urban surface is considered composed of:
		# 1 roof, 4 walls (2 of which are shadowed at any time) and 1 road
		# The typical urban unit is considered cubic 
		# so every length is multiplied by the roof width to get a surface
		# the urban fraction is used to derive the plan are of veg and built surfaces
		# the 3d components are then based on LAI and building height
	
		if(URBAN_FRAC>0) {
			TOTAL_URBAN_PLAN <- (ROAD_WIDTH+ROOF_WIDTH)*ROOF_WIDTH
			OVERALL_PLAN_AREA <- TOTAL_URBAN_PLAN/URBAN_FRAC
			TOTAL_URBAN_3D<-4*ROOF_HEIGHT*ROOF_WIDTH+ROAD_WIDTH*ROOF_WIDTH+ROOF_WIDTH*ROOF_WIDTH
			ACTIVE_WALL_3D<-2*active_wall*ROOF_HEIGHT*ROOF_WIDTH
			ACTIVE_ROOF_3D<-active_roof*ROOF_WIDTH*ROOF_WIDTH
			ACTIVE_ROAD_3D<-active_road*ROAD_WIDTH*ROOF_WIDTH
			TOTAL_ACTIVE_URBAN_3D<-ACTIVE_WALL_3D+ACTIVE_ROOF_3D+ACTIVE_ROAD_3D
		}else {
			OVERALL_PLAN_AREA <- 1
			TOTAL_URBAN_PLAN <- 0
			TOTAL_URBAN_3D <- 0
			ACTIVE_WALL_3D<- 0
			ACTIVE_ROOF_3D<- 0
			ACTIVE_ROAD_3D<- 0
			TOTAL_ACTIVE_URBAN_3D<- 0
		}
		#Vegetation fraction
		TOTAL_VEG_PLAN <- OVERALL_PLAN_AREA*(1-URBAN_FRAC)
		TOTAL_GRASS_PLAN <- TOTAL_VEG_PLAN*GRASS_FRAC
		TOTAL_TREE_PLAN <- TOTAL_VEG_PLAN*(1-GRASS_FRAC)

		ACTIVE_TREE_3D<-LAI_TREE*TOTAL_TREE_PLAN
		MAX_ACTIVE_TREE_3D<-LAI_MAX*TOTAL_TREE_PLAN
	
		# LAI constant for grass (people cut it!)
		ACTIVE_GRASS_3D<-LAI_GRASS*TOTAL_GRASS_PLAN
	
		TOTAL_ACTIVE_VEG_3D<-ACTIVE_TREE_3D+ACTIVE_GRASS_3D
		MAX_ACTIVE_VEG_3D<-MAX_ACTIVE_TREE_3D+max(LAI_GRASS,1)*TOTAL_GRASS_PLAN

		###############################################################
		# Compute needed indices
		# The idea here is to normalize by the total surface:
		# -> whole 3d urban surface + the max active vegetation surface
	      ###############################################################
	
		TOTAL_SURF_3D <- TOTAL_URBAN_3D+MAX_ACTIVE_VEG_3D
	
		ACTIVE_ROAD_INDEX <- ACTIVE_ROAD_3D/TOTAL_SURF_3D
		ACTIVE_WALL_INDEX <- ACTIVE_WALL_3D/TOTAL_SURF_3D
		ACTIVE_ROOF_INDEX <- ACTIVE_ROOF_3D/TOTAL_SURF_3D
		ACTIVE_URBAN_INDEX <- TOTAL_ACTIVE_URBAN_3D/TOTAL_SURF_3D
		#

		ACTIVE_VEG_INDEX <- array(data=TOTAL_ACTIVE_VEG_3D/TOTAL_SURF_3D,dim=length(lst_time))
		ACTIVE_TREE_INDEX <-array(data=ACTIVE_TREE_3D/TOTAL_SURF_3D,dim=length(lst_time))
		ACTIVE_GRASS_INDEX <- array(data=ACTIVE_GRASS_3D/TOTAL_SURF_3D,dim=length(lst_time))

		ACTIVE_TOTAL_INDEX <- (TOTAL_ACTIVE_URBAN_3D+TOTAL_ACTIVE_VEG_3D)/TOTAL_SURF_3D
		
	return(rbind(ACTIVE_URBAN_INDEX,ACTIVE_VEG_INDEX,ACTIVE_TOTAL_INDEX))

}

Output_Header <- function(FILE) {
	write("# FRAISE model output",file=FILE,append = TRUE, sep = " ")
	write("#",file=FILE,append = TRUE, sep = " ") 
	write("# SITE: a code for the location",file=FILE,append = TRUE, sep = " ")
	write("# DOY : day of year ",file=FILE,append = TRUE, sep = " ")
	write("# UZE : Urban Zone to characterize Energy partitioning (2 = Low Density; 3 = Medium Density; 4 = High Density)",file=FILE,append = TRUE, sep = " ")
	write("#",file=FILE,append = TRUE, sep = " ")
	write("# All following values are for mean day time (+/- 3 hours around noon)",file=FILE,append = TRUE, sep = " ")
	write("#",file=FILE,append = TRUE, sep = " ") 
	write("# CHI_TOT  : total active surface index                                                                    (-)",file=FILE,append = TRUE, sep = " ")
	write("# CHI_BUILT: active built index                                                                            (-)",file=FILE,append = TRUE, sep = " ")
	write("# CHI_VEG  : active vegetation index                                                                       (-)",file=FILE,append = TRUE, sep = " ")
	write("# QUP_RAT  : predicted ratio of outgoing radiant energy to that of incoming (excluding QF contribution)    (-)",file=FILE,append = TRUE, sep = " ")
	write("# QS_RAT   : predicted ratio of net storage heat (excluding QF contribution)                               (-)",file=FILE,append = TRUE, sep = " ")
	write("# QE_RAT   : predicted ratio of turbulent latent heat (excluding QF contribution)                          (-)",file=FILE,append = TRUE, sep = " ")
	write("# BOWEN    : predicted Bowen ratio                    (excluding QF contribution)                          (-)",file=FILE,append = TRUE, sep = " ")
	write("# QDOWN    : incoming raiant energy (user input)                                                          W m{-2}",file=FILE,append = TRUE, sep = " ")
	write("# QF       : anthropogenic heat (user input)                                                              W m{-2}",file=FILE,append = TRUE, sep = " ")
	write("# QUP      : predicted outgoing radiant flux (including QF contribution)                                  W m{-2}",file=FILE,append = TRUE, sep = " ")
	write("# QS       : predicted net storage heat flux (including QF contribution)                                  W m{-2}",file=FILE,append = TRUE, sep = " ")
	write("# QH_RES   : predicted (residual) turbulent sensible heat flux (including QF contribution)                W m{-2}",file=FILE,append = TRUE, sep = " ")
	write("# QH_BOW   : predicted (from Bowen ratio) turbulent sensible heat flux (including QF contribution)        W m{-2}",file=FILE,append = TRUE, sep = " ")
	write("# QE       : predicted turbulent latent heat flux (including QF contribution)                             W m{-2}",file=FILE,append = TRUE, sep = " ")
	write("# QSTAR    : predicted net all-wave radiative flux (including QF contribution)                            W m{-2}",file=FILE,append = TRUE, sep = " ")
	write("#",file=FILE,append = TRUE, sep = " ")
	write("# SITE    DOY   UZE  CHI_TOT  CHI_URB  CHI_VEG  QUP_RAT  QS_RAT   QE_RAT    BOWEN    QDOWN         QF        QUP        QS       QH_RES     QH_BOW        QE       QSTAR ",file=FILE,append = TRUE, sep = " ")
	write("#  -       -     -     -        -        -         -       -        -         -      W m{-2}     W m{-2}   W m{-2}    W m{-2}    W m{-2}    W m{-2}    W m{-2}    W m{-2}",file=FILE,append = TRUE, sep = " ")
	
}

