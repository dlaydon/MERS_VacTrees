require(plot3D)
require(here)
require(ggplot2)
#install.packages("hrbrthemes")
library(hrbrthemes)
library(grid)
library(gridExtra)
#install.packages("colorRamps")
library(colorRamps)
#install.packages("cowplot")
library("cowplot")

# Libraries
library(dplyr)
library(tidyr)
library(viridis)

ParDefaults = par()
OrigMAR 	= ParDefaults$mar

############################################################################
OrigMAI = par ("mai") ### good to record if using layout functions
OrigMAR = par ("mar") ### good to record if using layout functions

SummStat 				= function(Vector, IncludeMode = FALSE) ### returns summary statistics of vector
{
	Summary = c(mean(Vector), median(Vector), quantile(Vector, c(0.025,0.975), na.rm = T))
	NAMES 	= c("Mean", "Median", "LowerCrI", "UpperCrI")
	if (IncludeMode) 
	{
		Summary = c(Summary	, mlv(Vector, method = "mfv")$M)
		NAMES 	= c(NAMES	, "Mode")
	}
	names(Summary) = NAMES
	return (Summary)
}

ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R")
CppRootDirectory 	= file.path(ProjectDirectory, "MERS_Vac", "MERS_Vac") 
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
RawDataDirectory 	= file.path(ProjectDirectory, "Data"		)
source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))

# Import processed data for Cpp
DATA 		= read.table(file = file.path(RawDataDirectory, "MERS_forCpp.txt"), header = TRUE)
Day_0  		= as.Date("2013-01-01") #### 1st Jan 2013 is Day zero. 
Day_Final  	= Day_0 + max(DATA$onset) 

TotalDeaths 		= sum(DATA$Dead)
TotalDeaths_HCW 	= length(which(DATA$HCW == 1 & DATA$Dead == 1))
TotalDeaths_nHCW 	= length(which(DATA$HCW == 0 & DATA$Dead == 1))
LastDay_Data_Int 	= max(DATA$onset)

WindowsList = list(c(0, 548), c(0,180), c(181,364), c(365,548))

TotalDeaths_Window2 		= sum(DATA$Dead[which(							DATA$onset >= WindowsList[[2]][1] & DATA$onset <= WindowsList[[2]][2])])
TotalDeaths_HCW_Window2 	= length(which(DATA$HCW == 1 & DATA$Dead == 1 & DATA$onset >= WindowsList[[2]][1] & DATA$onset <= WindowsList[[2]][2]))
TotalDeaths_nHCW_Window2 	= length(which(DATA$HCW == 0 & DATA$Dead == 1 & DATA$onset >= WindowsList[[2]][1] & DATA$onset <= WindowsList[[2]][2]))

TotalDeaths_Window3 		= sum(DATA$Dead[which(							DATA$onset >= WindowsList[[3]][1] & DATA$onset <= WindowsList[[3]][2])])
TotalDeaths_HCW_Window3 	= length(which(DATA$HCW == 1 & DATA$Dead == 1 & DATA$onset >= WindowsList[[3]][1] & DATA$onset <= WindowsList[[3]][2]))
TotalDeaths_nHCW_Window3 	= length(which(DATA$HCW == 0 & DATA$Dead == 1 & DATA$onset >= WindowsList[[3]][1] & DATA$onset <= WindowsList[[3]][2]))

TotalDeaths_Window4 		= sum(DATA$Dead[which(							DATA$onset >= WindowsList[[4]][1] & DATA$onset <= WindowsList[[4]][2])])
TotalDeaths_HCW_Window4 	= length(which(DATA$HCW == 1 & DATA$Dead == 1 & DATA$onset >= WindowsList[[4]][1] & DATA$onset <= WindowsList[[4]][2]))
TotalDeaths_nHCW_Window4 	= length(which(DATA$HCW == 0 & DATA$Dead == 1 & DATA$onset >= WindowsList[[4]][1] & DATA$onset <= WindowsList[[4]][2]))

## make epidemic curve for Data
DATA_EpiCurves = rep(NA, LastDay_Data_Int + 1)
for (Day in 0:LastDay_Data_Int) DATA_EpiCurves[Day + 1] = length(which(DATA$onset == Day))
ConvertDailyToWeeklyInc = function(DailyInc)
{
	DaysToConsider 	= seq(1, length(DailyInc), by = 7)
	WeeklyInc 		= c(0, diff(cumsum(DailyInc), lag = 7))[DaysToConsider] # have to include zero for the first week 
	if (sum(WeeklyInc, na.rm = T) != sum(DailyInc, na.rm = T)) 
		WeeklyInc[length(WeeklyInc)] = sum(DailyInc, na.rm = T) - sum(WeeklyInc, na.rm = T)
	return(WeeklyInc)
}
DATA_EpiCurves_Weekly 			= ConvertDailyToWeeklyInc(DATA_EpiCurves)
DATA_EpidemicPeakSize_Weekly 	= max(DATA_EpiCurves_Weekly)

## Add #defines from Cpp code
AnimalReservoir			= -1 # note defined as 0 in Cpp but -1 here
SameHosptial			= 0 # i.e. transmission within same hospital
SameRegion				= 1 # i.e. transmission within same region but between different hospitals
DifferentRegion			= 2 # i.e. transmission between regions
AnimalReservoirLabelNum = 3

ChangeH2HLabel_AnimalTransmission = function(TreeData)
{
	TreeData$H2Htype[which(TreeData$H2Htype == DifferentRegion & TreeData$Infector == AnimalReservoir)] = AnimalReservoirLabelNum
	return(TreeData)
}
FixTrees = function(TreeData, BurninTrees = 900)
{
	TreeData 	= TreeData[TreeData$iteration >= BurninTrees, ]
	TreeData	= ChangeH2HLabel_AnimalTransmission(TreeData)
}

# contributions from all transmission types.
COLNAMES_TransContribMat_Props 	= c("Contrib_SameHosp", "Contrib_SameRegion", "Contrib_DiffRegion", "Contrib_Reservoir")
COLNAMES_TransContribMat_Abs 	= c("Contrib_SameHosp_MeanAbs", "Contrib_SameRegion_MeanAbs", "Contrib_DiffRegion_MeanAbs", "Contrib_Reservoir_MeanAbs")
COLNAMES_TransContribMat 		= c(COLNAMES_TransContribMat_Props, COLNAMES_TransContribMat_Abs)

Orig_CFs = c("PropCasesAverted", "PropCasesAverted_HCW", "PropCasesAverted_nHCW", 
		"DeathsAverted_CF", "DeathsAverted_CF_HCW", "DeathsAverted_CF_nHCW", "FinalCaseDate_CF", "Peak", "Peak_Reduction")
DirectIndirect_CFs = c("Cases_DirectlyProtected", "Cases_IndirectlyProtected", "Deaths_DirectlyProtected", "Deaths_IndirectlyProtected")
ListOfCounterfactuals = c(Orig_CFs, DirectIndirect_CFs)
StatNames = c("Mean", "LowerCrI", "UpperCrI")

## Make table of all model runs
ModelRuns = DefineModelRuns(
		
		IncludeCompletedRunsOnly 	= TRUE						, 
		
		VacCampStrategies 				= c("REACTIVE", "PROACTIVE"), 
#		ReactLevels						= c("HOSPITAL", "REGIONAL", "NATIONAL")							,
#		Efficacies_CamelControls		= seq(0, 0.5, 0.1)													,
		ImplementationDelays 			= seq(0, 28, 2)											, 
		
		WaningInReactiveAndNot 			= 1										, 
#		Pruning_Ranges 					= list(c(0, 548), c(0,180), c(181,364), c(365,548))			,
#		DoTriggersAndNot = 1	, 
#		TrigThresholds	= c(5,10,15,20,25,30)	,
#		TrigTimeframes	= c(7,14,30,60,90,120), 
#		ExpWaningAndNot					= 1:0														,
		
		ImmunityDelays 					= c(14, 0)													, 
		Efficacies_Start 				= seq(0.00, 1, by = 0.1)										, 
		VaccineDurations 				= c(0, 20, 15, 10, 5, 2, 1)														, 
		TimesSinceVaccination 			= c(0.5, 1:10)
) 
dim(ModelRuns)

TotalDeaths_LongVec			= rep(NA, dim(ModelRuns)[1])
TotalDeaths_HCW_LongVec		= rep(NA, dim(ModelRuns)[1])
TotalDeaths_nHCW_LongVec	= rep(NA, dim(ModelRuns)[1])

Window_1_Indices = which(ModelRuns$MinDay_Pruning == WindowsList[[1]][1] & ModelRuns$MaxDay_Pruning == WindowsList[[1]][2])
Window_2_Indices = which(ModelRuns$MinDay_Pruning == WindowsList[[2]][1] & ModelRuns$MaxDay_Pruning == WindowsList[[2]][2])
Window_3_Indices = which(ModelRuns$MinDay_Pruning == WindowsList[[3]][1] & ModelRuns$MaxDay_Pruning == WindowsList[[3]][2])
Window_4_Indices = which(ModelRuns$MinDay_Pruning == WindowsList[[4]][1] & ModelRuns$MaxDay_Pruning == WindowsList[[4]][2])

TotalDeaths_LongVec[Window_1_Indices] = TotalDeaths
TotalDeaths_LongVec[Window_2_Indices] = TotalDeaths_Window2
TotalDeaths_LongVec[Window_3_Indices] = TotalDeaths_Window3
TotalDeaths_LongVec[Window_4_Indices] = TotalDeaths_Window4

TotalDeaths_HCW_LongVec[Window_1_Indices] = TotalDeaths_HCW
TotalDeaths_HCW_LongVec[Window_2_Indices] = TotalDeaths_HCW_Window2
TotalDeaths_HCW_LongVec[Window_3_Indices] = TotalDeaths_HCW_Window3
TotalDeaths_HCW_LongVec[Window_4_Indices] = TotalDeaths_HCW_Window4

TotalDeaths_nHCW_LongVec[Window_1_Indices] = TotalDeaths_nHCW
TotalDeaths_nHCW_LongVec[Window_2_Indices] = TotalDeaths_nHCW_Window2
TotalDeaths_nHCW_LongVec[Window_3_Indices] = TotalDeaths_nHCW_Window3
TotalDeaths_nHCW_LongVec[Window_4_Indices] = TotalDeaths_nHCW_Window4

# remove ModelRuns
RedundantIndices = which(ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$ImmunityDelay == 0)
if (length(RedundantIndices) > 0)	ModelRuns = ModelRuns[-RedundantIndices, ]

dim(ModelRuns)
ModelRun = ModelRuns[1,]
ModelRuns$Pruning_Range = NULL

Category_Colnames = colnames(ModelRuns)  ## i.e. non-counterfactual cols
## Add counterfactuals to ModelRuns
for (counterfactualquantity in ListOfCounterfactuals)	
	for (Statistic in StatNames)
		ModelRuns[, paste0(counterfactualquantity, "_", Statistic)] = rep(NA, dim(ModelRuns)[1])

TransmissionContribMatrix 			= matrix(nrow = dim(ModelRuns), ncol = length(COLNAMES_TransContribMat))
colnames(TransmissionContribMatrix) = COLNAMES_TransContribMat

ModelRuns 				= cbind(ModelRuns, TransmissionContribMatrix)
Non_Category_Colnames 	= colnames(ModelRuns) [!(colnames(ModelRuns) %in% Category_Colnames)]

FilenamesCompleted 	= c()
DidRunFail 			= rep(0, dim(ModelRuns)[1])
RemoveBurnInPeriod 	= function(Chains, BurnIn = 1000) return (Chains[-(1:BurnIn), ])

dim(ModelRuns)
MR_index = 1
for (MR_index in 1:dim(ModelRuns)[1])
{
	ModelRun = ModelRuns[MR_index,]
	
	OutputString 	= paste0("_", ModelRun$OutputFolderNames)
	OutputSubDir 	= ModelRun$OutputFolderNames
	cat(paste0("MR ", MR_index, "/", dim(ModelRuns)[1]), " ", OutputSubDir, "\t")
	
	#### Import parameter and counterfactual chains. 
	CF_CHAINS_Summary = try(read.table(file = file.path(CppOutputDirectory, paste0("CF_ChainsSumm" 	, OutputString, ".txt")), header = T, sep = "\t"), silent = TRUE)
	if (class(CF_CHAINS_Summary) == "try-error")
	{
		CF_CHAINS 	= try(read.table(file = file.path(CppOutputDirectory, paste0("CF_Chains" 	, OutputString, ".txt")), header = T, sep = "\t"), silent = TRUE)
		NAMES 		= colnames(CF_CHAINS)
		
	} else NAMES = CF_CHAINS_Summary$Quantity
	if (class(CF_CHAINS_Summary) == "try-error" && class(CF_CHAINS) == "try-error") next
	
	if (exists("CF_CHAINS"))
	{
		#### Remove empty final columns
		while (all(is.na(CF_CHAINS	[, dim(CF_CHAINS)[2]]))) CF_CHAINS	[, dim(CF_CHAINS)[2]] = NULL
		#### Remove Burnin
		CF_CHAINS 		= RemoveBurnInPeriod(CF_CHAINS	, BurnIn = 1000)
	}
	
	DirIndirFileName = file.path(CppOutputDirectory, paste0("DirectIndirect_CF" 	, OutputString, ".txt"))
	DirectIndirect_Summary = try(read.table(file = DirIndirFileName, header = F, sep = "\t"), silent = TRUE)
	if (class(DirectIndirect_Summary) == "try-error") 
	{
		DidRunFail[MR_index] = 1
	
	} else {
		
		colnames(DirectIndirect_Summary) = StatNames
		rownames(DirectIndirect_Summary) = DirectIndirect_CFs
		
		counterfactualquantity = DirectIndirect_CFs[1]
		for (counterfactualquantity in DirectIndirect_CFs)	
		{
			Colnames_ThisCounterFactual 						= paste0(counterfactualquantity, "_", StatNames)
			ModelRuns[MR_index, Colnames_ThisCounterFactual] 	= DirectIndirect_Summary[counterfactualquantity, ]
		}
		
	}
	
	EpiCurvesFileName 		= file.path(CppOutputDirectory, paste0("CF_EpiCurves" 			, OutputString, ".txt"))
	TreesFileName			= file.path(CppOutputDirectory, paste0("Trees_CF"				, OutputString, ".txt"))
	
	if (file.exists(EpiCurvesFileName))
	{
		CF_EpiCurves 	= try(read.table(file = EpiCurvesFileName, header = T, sep = "\t"), silent = TRUE)
		if (class(CF_EpiCurves) != "try-error")
			if (dim(CF_EpiCurves)[1] > 100) 
			{
				cat (paste0("MR_index ", MR_index, " ", OutputString, " EpiCurves not weekly"))
				DidRunFail[MR_index] = 1
			}
	} else 
	{
		cat (paste0("MR_index ", MR_index, " ", OutputString, " EpiCurves don't exist"))
		DidRunFail[MR_index] = 1
	}
	
	
	counterfactualquantity = Orig_CFs[1]
	for (counterfactualquantity in Orig_CFs)	
	{
		if (counterfactualquantity %in% NAMES)
		{
			Colnames_ThisCounterFactual 						= paste0(counterfactualquantity, "_", StatNames)
			if (exists("CF_CHAINS"))
			{
				Summ_CounterFactual									= signif(SummStat(CF_CHAINS[, counterfactualquantity]), 3)
				ModelRuns[MR_index, Colnames_ThisCounterFactual] 	= Summ_CounterFactual[StatNames]
				rm(Summ_CounterFactual)
				
			} else ModelRuns[MR_index, Colnames_ThisCounterFactual] = CF_CHAINS_Summary[CF_CHAINS_Summary$Quantity == counterfactualquantity, StatNames]
			
		} else if (counterfactualquantity == "Peak")
		{
			if (file.exists(EpiCurvesFileName))
			{
				if (class(CF_EpiCurves) != "try-error")
				{
					### Convert
					for (Statistic in StatNames)
					{
						Colname_ThisCounterFactual 						= paste0(counterfactualquantity, "_", Statistic)
						Peak_ThisStatistic 								= max(CF_EpiCurves[, Statistic])
						ModelRuns[MR_index, Colname_ThisCounterFactual] = Peak_ThisStatistic
						
					}					
				}
				
			} 
		}
	}
	
	if (file.exists(TreesFileName))
	{
		CF_Trees = try(read.table(file = TreesFileName, header = T), silent = TRUE)
		if (class(CF_Trees) != "try-error")
		{
			# process tree (take out burnin and differentiate properly between animal reservoir and between region transmission)
			CF_Trees = FixTrees(CF_Trees)
			
			if (dim(CF_Trees)[1] > 0)
			{
				TableDummy = table(CF_Trees$H2Htype) 
				if (length(TableDummy) == 4)
				{
					ModelRuns[MR_index, COLNAMES_TransContribMat_Props] = signif(TableDummy / dim(CF_Trees)[1], 3)
					ModelRuns[MR_index, COLNAMES_TransContribMat_Abs] 	= round(colSums(table(CF_Trees$iteration, CF_Trees$H2Htype)) / length(unique(CF_Trees$iteration)))
					
				} else	DidRunFail[MR_index] = 1
				
			} else 
			{
				cat (paste0("\nMR_index ", MR_index, " ", OutputString, " CF_Trees have too few iterations"))
				DidRunFail[MR_index] = 1
			}
		}
		
	} else 
	{
		cat (paste0("MR_index ", MR_index, " ", OutputString, " CF_Trees don't exist"))
		DidRunFail[MR_index] = 1
	}
	
	# record ModelRun again after having added counterfactual information.
	ModelRun = ModelRuns[MR_index,]
	
	#### clean up. 
	FilenamesCompleted = c(FilenamesCompleted, OutputString)
	rm(CF_CHAINS, CF_CHAINS_Summary, OutputSubDir, CF_EpiCurves, CF_EpiCurves_Deaths)
	cat("\n")
}

# Post process deaths because in your infinite wisdom you weren't consistent between cases and deaths in Cpp code.
for (Statistic in StatNames)
{
	ModelRuns[, paste0("DeathsAverted_CF"		, "_", Statistic)] = ModelRuns[, paste0("DeathsAverted_CF"		, "_", Statistic)] / TotalDeaths
	ModelRuns[, paste0("DeathsAverted_CF_HCW"	, "_", Statistic)] = ModelRuns[, paste0("DeathsAverted_CF_HCW"	, "_", Statistic)] / TotalDeaths_HCW
	ModelRuns[, paste0("DeathsAverted_CF_nHCW"	, "_", Statistic)] = ModelRuns[, paste0("DeathsAverted_CF_nHCW"	, "_", Statistic)] / TotalDeaths_nHCW
}
ModelRuns$VaccineDuration[which(ModelRuns$VaccineDuration == 0)] = Inf

## Convert "Peak" counterfactuals to "reduction in peak" counterfactuals.
ModelRuns$Peak_Reduction_Mean 		= (DATA_EpidemicPeakSize_Weekly - ModelRuns$Peak_Mean		) 	/ DATA_EpidemicPeakSize_Weekly
ModelRuns$Peak_Reduction_LowerCrI 	= (DATA_EpidemicPeakSize_Weekly - ModelRuns$Peak_UpperCrI	)	/ DATA_EpidemicPeakSize_Weekly # note switch of upper and lower CrI
ModelRuns$Peak_Reduction_UpperCrI 	= (DATA_EpidemicPeakSize_Weekly - ModelRuns$Peak_LowerCrI	)	/ DATA_EpidemicPeakSize_Weekly # note switch of upper and lower CrI

# write.table
write.table(ModelRuns, file = file.path(ProjectDirectory, paste0("ModelRunsSummary.txt")), row.names = F, col.names = T, quote = F, sep = "\t")




