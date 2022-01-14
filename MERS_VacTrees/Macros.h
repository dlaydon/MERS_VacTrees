#pragma once

//#define PRINT_PROGRAM_PROGRESS
#define GET_VARIABLE_NAME(Variable) (#Variable)

#define AnimalReservoir		0

#define Alive				0
#define Dead				1

#define Asymptomatic		0
#define Symptomatic			1

//// H2Htype
#define SameHosptial		0
#define SameRegion			1
#define DifferentRegion		2

//// _spatialLevel
#define SL_No_Space				0 // _spatialLevel 0: no space / 1:region / 2: hospital/ 3: region & hospital
#define SL_Region				1 // _spatialLevel 0: no space / 1:region / 2: hospital/ 3: region & hospital
#define	SL_Hospital				2 // _spatialLevel 0: no space / 1:region / 2: hospital/ 3: region & hospital
#define SL_RegionAndHospital	3 // _spatialLevel 0: no space / 1:region / 2: hospital/ 3: region & hospital

/// _variabilityRAtClusterLevel
#define vC_No_Hetero				0	//For _variabilityRAtClusterLevel - 0: no heterogeneity / 1: individual level / 2: cluster level
#define vC_IndividualLevel			1	//For _variabilityRAtClusterLevel - 0: no heterogeneity / 1: individual level / 2: cluster level
#define vC_ClusterLevel				2	//For _variabilityRAtClusterLevel - 0: no heterogeneity / 1: individual level / 2: cluster level

#define No_Seasonality				0	// 0: no seasonality	; 1: seasonality (effect seasonality can't be below 1)	; 2: true cos function
#define Seasonality_AboveOne		1	// 0: no seasonality	; 1: seasonality (effect seasonality can't be below 1)	; 2: true cos function
#define Seasonality_TrueCos			2	// 0: no seasonality	; 1: seasonality (effect seasonality can't be below 1)	; 2: true cos function

#define MD_Value					9999
#define MAX_ONSET_DAY				549		

#define MEAN		0
#define MEDIAN		1
#define LOWER_CrI	2
#define UPPER_CrI	3