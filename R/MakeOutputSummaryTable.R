rm(list=ls())
Sys.time()
gc()

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
options(width = 108L)
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
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppRootDirectory 	= file.path(ProjectDirectory, "MERS_Vac"	, "MERS_Vac") 
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
RawDataDirectory 	= file.path(ProjectDirectory, "Data"		)
source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))

# Import processed data for Cpp
DATA 		= read.table(file = file.path(RawDataDirectory, "DJL_MERS_forCpp.txt"), header = TRUE)
Day_0  		= as.Date		("2013-01-01") #### 1st Jan 2013 is Day zero. 
Day_Final  	= Day_0 + max(DATA$onset) # note Different from Day_Final in MakeLineListForCpp.R script.

TotalDeaths 		= sum(DATA$Dead)
TotalDeaths_HCW 	= length(which(DATA$HCW == 1 & DATA$Dead == 1))
TotalDeaths_nHCW 	= length(which(DATA$HCW == 0 & DATA$Dead == 1))
LastDay_Data_Int 	= max(DATA$onset)

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
DATA_EpiCurves_Weekly = ConvertDailyToWeeklyInc(DATA_EpiCurves)
DATA_EpidemicPeakSize_Weekly = max(DATA_EpiCurves_Weekly)

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
ListOfCounterfactuals = c("PropCasesAverted", "PropCasesAverted_HCW", "PropCasesAverted_nHCW", 
		"DeathsAverted_CF", "DeathsAverted_CF_HCW", "DeathsAverted_CF_nHCW", "FinalCaseDate_CF", "Peak", "Peak_Reduction")
StatNames = c("Mean", "LowerCrI", "UpperCrI")

## Make table of all model runs
ModelRuns = DefineModelRuns(
		
		IncludeCompletedRunsOnly = TRUE, 
		
		VacCampStrategies 				= c("REACTIVE", "PROACTIVE"), 
		ReactLevels						= c("HOSPITAL", "REGIONAL", "NATIONAL")							,
		Efficacies_CamelControls		= seq(0, 0.5, 0.1)													,
		ImplementationDelays 			= seq(0, 28, 2)											, 
		ImmunityDelays 					= c(14, 0)													, 
		Efficacies_Start 				= seq(0.00, 1, by = 0.05)										, 
		VaccineDurations 				= c(0, 20, 15, 10, 5, 2, 1)														, 
		TimesSinceVaccination 			= c(0.5, 1:10)
) 

# remove archaic ModelRuns
RedundantIndices = which(ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$ImmunityDelay == 0)
if (length(RedundantIndices) > 0)	ModelRuns = ModelRuns[-RedundantIndices, ]

dim(ModelRuns)
ModelRuns$OutputFolderNames
ModelRun = ModelRuns[1,]

Category_Colnames = colnames(ModelRuns)  ## i.e. non-counterfactual cols
## Add counterfactuals to ModelRuns
for (counterfactualquantity in ListOfCounterfactuals)	
	for (Statistic in StatNames)
		ModelRuns[, paste0(counterfactualquantity, "_", Statistic)] = rep(NA, dim(ModelRuns)[1])

TransmissionContribMatrix 			= matrix(nrow = dim(ModelRuns), ncol = length(COLNAMES_TransContribMat))
colnames(TransmissionContribMatrix) = COLNAMES_TransContribMat

FilenamesCompleted = c()
DidRunFail 			= rep(0, dim(ModelRuns)[1])
RemoveBurnInPeriod 	= function(Chains, BurnIn = 1000) return (Chains[-(1:BurnIn), ])


for (MR_index in 1:dim(ModelRuns)[1])
{
	ModelRun = ModelRuns[MR_index,]
	
	OutputString 	= ChooseOutputString(ModelRun, Folder = FALSE	)
	OutputSubDir 	= ChooseOutputString(ModelRun, Folder = TRUE	)
	cat(paste0("MR ", MR_index, "/", dim(ModelRuns)[1]), " ", OutputSubDir, "\t")
	
	#### Import parameter and counterfactual chains. 
	CF_CHAINS 	= read.table(file = file.path(CppOutputDirectory, paste0("CF_Chains" 	, OutputString, ".txt")), header = T, sep = "\t")	
	#### Remove empty final columns
	while (all(is.na(CF_CHAINS	[, dim(CF_CHAINS)[2]]))) CF_CHAINS	[, dim(CF_CHAINS)[2]] = NULL
	
	#### Remove Burnin
	CF_CHAINS 		= RemoveBurnInPeriod(CF_CHAINS	, BurnIn = 1000)
	
	EpiCurvesFileName 		= file.path(CppOutputDirectory, paste0("CF_EpiCurves" 			, OutputString, ".txt"))
	EpiCurvesDeathsFileName = file.path(CppOutputDirectory, paste0("CF_EpiCurves_Deaths" 	, OutputString, ".txt"))
	TreesFileName			= file.path(CppOutputDirectory, paste0("Trees_CF"				, OutputString, ".txt"))
	
	if (file.exists(EpiCurvesFileName))
	{
		CF_EpiCurves 	= read.table(file = EpiCurvesFileName, header = T, sep = "\t")
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
	
	if (file.exists(EpiCurvesDeathsFileName))
	{
		CF_EpiCurves_Deaths 	= read.table(file = EpiCurvesDeathsFileName, header = T, sep = "\t")
		if (dim(CF_EpiCurves_Deaths)[1] > 100) 
		{
			cat (paste0("MR_index ", MR_index, " ", OutputString, " CF_EpiCurves_Deaths not weekly"))
			DidRunFail[MR_index] = 1
		}
		
	} else 
	{
		cat (paste0("MR_index ", MR_index, " ", OutputString, " EpiCurves_Deaths don't exist"))
		DidRunFail[MR_index] = 1
	}
	
	colnames(CF_CHAINS)
	
	
	counterfactualquantity = ListOfCounterfactuals[1]
	for (counterfactualquantity in ListOfCounterfactuals)	
	{
		if (counterfactualquantity %in% colnames(CF_CHAINS))
		{
			Summ_CounterFactual	= signif(SummStat(CF_CHAINS[, counterfactualquantity]), 3)
			Colnames_ThisCounterFactual = paste0(counterfactualquantity, "_", StatNames)
			ModelRuns[MR_index, Colnames_ThisCounterFactual] = Summ_CounterFactual[StatNames]
			rm(Summ_CounterFactual)
			
		} else if (counterfactualquantity == "Peak")
		{
			if (file.exists(EpiCurvesFileName))
			{
				### Convert
				for (Statistic in StatNames)
				{
					Colname_ThisCounterFactual 	= paste0(counterfactualquantity, "_", Statistic)
					Peak_ThisStatistic = max(CF_EpiCurves[, Statistic])
					ModelRuns[MR_index, Colname_ThisCounterFactual] = Peak_ThisStatistic
					
				}
			} else cat (paste0("MR_index ", MR_index, " ", OutputString, " EpiCurves don't exist"))
		}
	}
	
	if (file.exists(TreesFileName))
	{
		CF_Trees = read.table(file = TreesFileName, header = T)
		# process tree (take out burnin and differentiate properly between animal reservoir and between region transmission)
		CF_Trees = FixTrees(CF_Trees)
		
		TransmissionContribMatrix[MR_index, COLNAMES_TransContribMat_Props] = signif(table(CF_Trees$H2Htype) / dim(CF_Trees)[1], 3)
		TransmissionContribMatrix[MR_index, COLNAMES_TransContribMat_Abs] 	= round(colSums(table(CF_Trees$iteration, CF_Trees$H2Htype)) / length(unique(CF_Trees$iteration)))
		
	} else 
	{
		cat (paste0("MR_index ", MR_index, " ", OutputString, " CF_Trees don't exist"))
		DidRunFail[MR_index] = 1
	}
	
	# record ModelRun again after having added counterfactual information.
	ModelRun = ModelRuns[MR_index,]
	
	#### clean up. 
	FilenamesCompleted = c(FilenamesCompleted, OutputString)
	rm(CF_CHAINS, OutputSubDir, CF_EpiCurves, CF_EpiCurves_Deaths)
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

ModelRuns = cbind(ModelRuns, TransmissionContribMatrix)

# write.table
write.table(ModelRuns, file = file.path(CppOutputDirectory, paste0("ModelRunsSummary.txt")), row.names = F, col.names = T, quote = F, sep = "\t")




