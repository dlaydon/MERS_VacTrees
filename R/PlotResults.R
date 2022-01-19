### Superceded by various other smaller leaner scripts - use these from now on.
##### Plots results from various model runs. 
##### Code nicked from ProducePlots9.R script in Dengue project. 
##### However structure different in important ways. Instead of "taxonomy" of files being "Kind of Plot" (e.g. posteriors) then "Model Run", will be the other way around to make model comparison easier.

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


IMPORT_SUMMARY_TABLE = TRUE
#IMPORT_SUMMARY_TABLE = FALSE

#MAKE_PLOTS = TRUE
MAKE_PLOTS = FALSE

MAKE_SUMMARY_TABLE = TRUE
#MAKE_SUMMARY_TABLE = FALSE

WEEKLY_INCIDENCE = TRUE
#WEEKLY_INCIDENCE = FALSE

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
Make95CrI_Char			= function(Mean, Lower, Upper)
{
	return(as.factor(paste0(Mean, " (", Lower, ", ", Upper, ")")))
}		
Make95CrI_CharVecVerison			= function(VectorOf_Mean_Lower_Upper)
{
	return(as.factor(paste0(VectorOf_Mean_Lower_Upper["Mean"], " (", VectorOf_Mean_Lower_Upper["LowerCRI"], ", ", VectorOf_Mean_Lower_Upper["UpperCrI"], ")")))
}

PosteriorCol 	= "turquoise"
PNG_res 		= 300

ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppRootDirectory 	= file.path(ProjectDirectory, "MERS_Vac"	, "MERS_Vac") 
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
RawDataDirectory 	= file.path(ProjectDirectory, "Data"		)
PlotsDirectory		= file.path(ProjectDirectory, "Plots"		)

source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))
source(file.path(R_ScriptDirectory, "LayoutPlottingFunctions.R"))

# Import processed data for Cpp
DATA 		= read.table(file = file.path(RawDataDirectory, "DJL_MERS_forCpp.txt"), header = TRUE)
Day_0  		= as.Date		("2013-01-01") #### 1st Jan 2013 is Day zero. 
Day_Final  	= Day_0 + max(DATA$onset) # note Different from Day_Final in MakeLineListForCpp.R script.
LastDay_Data_Int 	= max(DATA$onset)

TotalDeaths 		= sum(DATA$Dead)
TotalDeaths_HCW 	= length(which(DATA$HCW == 1 & DATA$Dead == 1))
TotalDeaths_nHCW 	= length(which(DATA$HCW == 0 & DATA$Dead == 1))

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

# import who-infected-whom trees without counterfactuals
Trees_noCFs = read.table(here("Output/ExampleTrees_nonCF.txt"), header = TRUE)

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
BurninTrees = 900
diff(range(Trees_noCFs$iteration))
Trees_noCFs	= FixTrees(Trees_noCFs)
diff(range(Trees_noCFs$iteration))
# contributions from all transmission types.
signif(table(Trees_noCFs$H2Htype) / dim(Trees_noCFs)[1], 3)
sum(table(Trees_noCFs$H2Htype) / dim(Trees_noCFs)[1] ) == 1



#COLNAMES_TransContribMat = c("Contrib_SameHosp", "Contrib_SameRegion", "Contrib_DiffRegion", "Contrib_Reservoir")

COLNAMES_TransContribMat_Props 	= c("Contrib_SameHosp", "Contrib_SameRegion", "Contrib_DiffRegion", "Contrib_Reservoir")
COLNAMES_TransContribMat_Abs 	= c("Contrib_SameHosp_MeanAbs", "Contrib_SameRegion_MeanAbs", "Contrib_DiffRegion_MeanAbs", "Contrib_Reservoir_MeanAbs")
COLNAMES_TransContribMat 		= c(COLNAMES_TransContribMat_Props, COLNAMES_TransContribMat_Abs)
ListOfCounterfactuals = c("PropCasesAverted", "PropCasesAverted_HCW", "PropCasesAverted_nHCW", 
		"DeathsAverted_CF", "DeathsAverted_CF_HCW", "DeathsAverted_CF_nHCW", "FinalCaseDate_CF", "Peak", "Peak_Reduction")
StatNames = c("Mean", "LowerCrI", "UpperCrI")

CF_Colnames = paste0(rep(ListOfCounterfactuals, each = 3), "_", StatNames)


if (IMPORT_SUMMARY_TABLE & !MAKE_PLOTS)
{
	ModelRuns = read.table(file = file.path(CppOutputDirectory, paste0("ModelRunsSummary.txt")), header = T, sep = "\t") 
				
} else {
	
	ModelRuns = DefineModelRuns(
			
			IncludeCompletedRunsOnly = TRUE, 
			
			VacCampStrategies 				= c("REACTIVE", "PROACTIVE"), 
			#VacCampStrategies 				= c("REACTIVE"), 
			#VacCampStrategies 				= c("PROACTIVE"), 
			ReactLevels						= c("HOSPITAL", "REGIONAL", "NATIONAL")							,
			#ReactLevels						= c("HOSPITAL")											,
			Efficacies_CamelControls		= seq(0, 0.5, 0.1)													,
			ImplementationDelays 			= seq(0, 28, 2)											, 
			ImmunityDelays 					= c(14, 0)													, 
			Efficacies_Start 				= seq(0.00, 1, by = 0.05)										, 
			VaccineDurations 				= c(0, 20, 15, 10, 5, 2, 1)														, 
			TimesSinceVaccination 			= c(0.5, 1:10)
	) 
	
	# remove hangover ModelRuns that are archaic and mess everything up. 
	RedundantIndices = which(ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$ImmunityDelay == 0)
	if (length(RedundantIndices) > 0)	ModelRuns 	= ModelRuns[-RedundantIndices, ]
	
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
	PlotChain 			= function(Chains = CHAINS, ParamNum = 1, N_Signif = 3, ...)
	{
		if (ParamNum <= dim(Chains)[2]) ## else do nothing. 
		{
			ChainToPlot = Chains[, ParamNum]
			ParamString	= paste0(ParamNum, ": ",colnames(Chains)[ParamNum])
			SummParam	= signif(SummStat(ChainToPlot), N_Signif)
			
			PlotTitle 	= paste0(ParamString, "\n", SummParam["Mean"], " (", SummParam["LowerCrI"]," - ", SummParam["UpperCrI"], ")")
			
			plot(1:length(ChainToPlot), ChainToPlot, type = "l", ylab = "", xlab = "Iteration", main = PlotTitle,	...)
		}
	}
	PlotPosterior 		= function(Chains = CHAINS, ParamNum = 1, N_Signif = 3, ...)
	{	
		if (ParamNum <= dim(Chains)[2]) ## else do nothing. 
		{
			ChainToPlot = Chains[, ParamNum]
			ParamString	= paste0(ParamNum, ": ", colnames(Chains)[ParamNum])
			SummParam	= signif(SummStat(ChainToPlot), N_Signif)
			PlotTitle 	= paste0(ParamString, "\n", SummParam["Mean"], " (", SummParam["LowerCrI"]," - ", SummParam["UpperCrI"], ")")
			hist(ChainToPlot, xlab = colnames(Chains)[ParamNum], freq = FALSE, col = PosteriorCol, main = PlotTitle, ...)
		}
	}
	
	RemoveBurnInPeriod 	= function(Chains, BurnIn = 1000) return (Chains[-(1:BurnIn), ])
	
	#head(ModelRuns[, Category_Colnames])
	DidRunFail 	= rep(0, dim(ModelRuns)[1])
	MR_index 	= 62
	
	
	for (MR_index in 1:dim(ModelRuns)[1])
	#for (MR_index in MR_index:dim(ModelRuns)[1])
	#for (MR_index in which(DidRunFail == 1))
	{
		## close any open plot devices
		CloseOpenPlotDevices()	
		ModelRun = ModelRuns[MR_index,]
		
		OutputString 	= ChooseOutputString(ModelRun, Folder = FALSE	)
		OutputSubDir 	= ChooseOutputString(ModelRun, Folder = TRUE	)
		cat(paste0("MR ", MR_index, "/", dim(ModelRuns)[1]), " ", OutputSubDir, "\t")
		
		#### Import parameter and counterfactual chains. 
		#CHAINS 		= read.table(file = file.path(CppOutputDirectory, paste0("Chains" 		, OutputString, ".txt")), header = T, sep = "\t", check.names = FALSE) ## check.names = FALSE so that underscore prefixes not replaced with X_	 
		CF_CHAINS 	= read.table(file = file.path(CppOutputDirectory, paste0("CF_Chains" 	, OutputString, ".txt")), header = T, sep = "\t")	
		#### Remove empty final columns
		#while (all(is.na(CHAINS		[, dim(CHAINS	)[2]]))) CHAINS		[, dim(CHAINS	)[2]] = NULL
		while (all(is.na(CF_CHAINS	[, dim(CF_CHAINS)[2]]))) CF_CHAINS	[, dim(CF_CHAINS)[2]] = NULL
		
		#### Remove Burnin
		#CHAINS 		= RemoveBurnInPeriod(CHAINS		, BurnIn = 1000)
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
			## check 
			#all(CF_EpiCurves[, "UpperCrI"] >= CF_EpiCurves[, "LowerCrI"])
			#all(CF_EpiCurves[, "UpperCrI"] >= CF_EpiCurves[, "Mean"])
			#all(CF_EpiCurves[, "Mean"] >= CF_EpiCurves[, "LowerCrI"])
			#
			#cbind(CF_EpiCurves[, "Mean"], CF_EpiCurves[, "LowerCrI"], CF_EpiCurves[, "Mean"] >= CF_EpiCurves[, "LowerCrI"])
			
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
			## check 
			#all(CF_EpiCurves_Deaths[, "UpperCrI"] >= CF_EpiCurves[, "LowerCrI"])
			#all(CF_EpiCurves_Deaths[, "UpperCrI"] >= CF_EpiCurves[, "Mean"])
			#all(CF_EpiCurves_Deaths[, "Mean"] >= CF_EpiCurves[, "LowerCrI"])
			#
			#cbind(CF_EpiCurves_Deaths[, "Mean"], CF_EpiCurves_Deaths[, "LowerCrI"], CF_EpiCurves_Deaths[, "Mean"] >= CF_EpiCurves_Deaths[, "LowerCrI"])
			
		} else 
		{
			cat (paste0("MR_index ", MR_index, " ", OutputString, " EpiCurves_Deaths don't exist"))
			DidRunFail[MR_index] = 1
		}
		
		colnames(CF_CHAINS)
		
		if (MAKE_SUMMARY_TABLE)
		{
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
						#plot(DATA_EpiCurves_Weekly, col = "black", ylab = "cases", xlab = "week", type = "l", lwd = 3)
						### Convert
						for (Statistic in StatNames)
						{
							#WeeklyInc 					= ConvertDailyToWeeklyInc(CF_EpiCurves[, Statistic])
							Colname_ThisCounterFactual 	= paste0(counterfactualquantity, "_", Statistic)
							Peak_ThisStatistic = max(CF_EpiCurves[, Statistic])
							ModelRuns[MR_index, Colname_ThisCounterFactual] = Peak_ThisStatistic
							#if (Statistic == "Mean") COL = "red" else COL = "pink"
							#points(WeeklyInc, col = COL, type = "l", lwd = 3)
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
		}
		
		# record ModelRun again after having added counterfactual information.
		ModelRun = ModelRuns[MR_index,]
		
		if (MAKE_PLOTS)
		{
			ModelRunPlot_dir = file.path(PlotsDirectory, OutputSubDir) ### given only one file, perhaps better to have all the files together, not subdivided by folder. 
			dir.create(ModelRunPlot_dir, showWarnings = FALSE)
			
			# later add this to posteriors
			if (ModelRun$VacCampStrategy == "PROACTIVE")
			{
				png(file = file.path(ModelRunPlot_dir, "EfficacyOverTime.png"), res = PNG_res, units = "in", width = 5, height = 5) ## width and height chosen to be dimensions of one screen. 
				YearsPostDose = 0:3650 / 365
				plot (YearsPostDose, ModelRun$Efficacy_Start * exp(-YearsPostDose / ModelRun$VaccineDuration), ylim = c(0, 1), 
						ylab = "Efficacy", xlab = "Time (years)", 
						main = paste0("Effiacy over time: Initial = ", ModelRun$Efficacy_Start*100, "%\nmean duration = ", ModelRun$VaccineDuration, " years"),  
						col = "turquoise")
				dev.off()
			}
			
			#### Plot all CF_chains
			cat(paste0("CF_Chains\t"))
			NumRows = 2; NumCols = 3; NumChainPlotsPerPNG = NumRows * NumCols; 
			png(file = file.path(ModelRunPlot_dir, "CF_ParamChains.png"), res = PNG_res, units = "in", width = 21, height = 12) ## width and height chosen to be dimensions of one screen. 
			SetUpMultiPlot(NoCols = NumCols, NoRows = NumRows, MultiPlot_Title = paste0("CF_ParamChains: ", OutputSubDir), Ratio_Plots_To_Title = 6, BY_ROW = TRUE)
			for (ParamNum in 1:NumChainPlotsPerPNG) 
				PlotChain(Chains = CF_CHAINS, ParamNum = ParamNum, cex.main = 2, cex.axis = 2, cex.lab = 2)
			dev.off()
			
			#### Plot all CF_posteriors
			cat(paste0("CF_Posteriors\t"))
			NumRows = 2; NumCols = 3; NumChainPlotsPerPNG = NumRows * NumCols; 
			png(file = file.path(ModelRunPlot_dir, "CF_ParamPosts.png"), res = PNG_res, units = "in", width = 21, height = 12) ## width and height chosen to be dimensions of one screen. 
			SetUpMultiPlot(NoCols = NumCols, NoRows = NumRows, 
					MultiPlot_Title = paste0("CF_ParamChains: ", OutputSubDir), 
					Ratio_Plots_To_Title = 6, BY_ROW = TRUE)
			for (ParamNum in 1:NumChainPlotsPerPNG) 
				PlotPosterior(Chains = CF_CHAINS, ParamNum = ParamNum, cex.main = 2, cex.axis = 2, cex.lab = 2)
			dev.off()
			
			
			### counterfactual epidemic curves
			DataCol = "blue"
			CF_Col = "red"; CF_Col_CrIs = "pink"
			png(file = file.path(ModelRunPlot_dir, "CF_EpiCurves.png"), res = PNG_res, width = 5.5, height = 5, units = "in")
		
			Xaxis = Day_0 + 0:max(DATA$onset)
			plot(Xaxis, DATA_EpiCurves, type = "h", xaxt = "n", col = DataCol, lwd = 3,
					xlab = "", ylab = "Cases",
					xaxt = "n", 
					xaxp  = c(0, 110, 10), 
					main = paste0("Daily MERS-nCoV cases in 2013-2014 KSA outbreak\n", OutputSubDir))
			points(Xaxis, CF_EpiCurves$UpperCrI, col = CF_Col_CrIs, xaxt = "n", lwd = 3, type = "h")
			points(Xaxis, CF_EpiCurves$Mean, col = CF_Col, xaxt = "n", lwd = 3, type = "h")
			points(Xaxis, CF_EpiCurves$LowerCrI, col = CF_Col_CrIs, xaxt = "n", lwd = 3, type = "h")
			points(Xaxis, DATA_EpiCurves, col = DataCol, xaxt = "n", lwd = 0.2)
			axis.Date(1, at = seq(Day_0, Day_Final, by = "months"), format = "%b-%y", las = 2, srt = 35)
			legend("topleft", legend = c("data", "counterfactual mean", "counterfactual CrIs"), col = c(DataCol, CF_Col, CF_Col_CrIs), pch = c(1,15, 15), pt.cex = 2)
			
			dev.off()
		}
		
		#### clean up. 
		FilenamesCompleted = c(FilenamesCompleted, OutputString)
		rm(CHAINS, CF_CHAINS, OutputSubDir, NumRows, NumCols, NumChainPlotsPerPNG, CF_EpiCurves, CF_EpiCurves_Deaths)
		CloseOpenPlotDevices()	
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
	
	DuffRuns = which(is.na(ModelRuns$PropCasesAverted_nHCW_Mean) | is.na(ModelRuns$Peak_Mean))
	which(DidRunFail == 1)
	sum(DidRunFail)
	length(DuffRuns)
	if (length(DuffRuns) != 0)
	{
		ModelRuns[DuffRuns, Category_Colnames] # print only categories.
		# save duffruns and run them again. 
		write.table(ModelRuns[DuffRuns, Category_Colnames], file = file.path(CppOutputDirectory, paste0("MissedRuns_291110.txt")), row.names = F, col.names = T, quote = F, sep = "\t")
	}
	summary(ModelRuns[DuffRuns, Category_Colnames] )
	
	# write.table
	write.table(ModelRuns, file = file.path(CppOutputDirectory, paste0("ModelRunsSummary.txt")), row.names = F, col.names = T, quote = F, sep = "\t")
	
}
	
head(ModelRuns[,CF_Colnames])
colnames(ModelRuns)
Category_Colnames = colnames(ModelRuns)[!(colnames(ModelRuns) %in% CF_Colnames) & !(colnames(ModelRuns) %in% COLNAMES_TransContribMat)]  ## i.e. non-counterfactual cols


# Print various quantities that you refer to in the text to the console

PrintCols = c("PropCasesAverted_Mean", "PropCasesAverted_LowerCrI", "PropCasesAverted_UpperCrI", 
		"DeathsAverted_CF_Mean", "DeathsAverted_CF_LowerCrI", "DeathsAverted_CF_UpperCrI")

ModelRuns[ModelRuns$VacCampStrategy == "PROACTIVE" & 
				ModelRuns$Efficacy_Start == 0.9 & 
				ModelRuns$VaccineDuration == 20 & 
				ModelRuns$TimeSinceVaccination %in% c(0.5) & 
				ModelRuns$Efficacy_CamelControls == unique(ModelRuns$Efficacy_CamelControls)[6], 	PrintCols]

unique(ModelRuns$Efficacy_CamelControls) == 0.30
class(ModelRuns$Efficacy_CamelControls)

ModelRuns[ModelRuns$VacCampStrategy == "PROACTIVE" & 
				ModelRuns$Efficacy_Start == 0.25 & 
				ModelRuns$VaccineDuration == Inf & 
				ModelRuns$TimeSinceVaccination %in% 0.5 & 
				ModelRuns$Efficacy_CamelControls == 0.0			, 	PrintCols]



ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & 
				ModelRuns$ReactLevel == "HOSPITAL" &
				ModelRuns$Efficacy_Start == 0.5 & 
				ModelRuns$ImplementationDelay %in% c(0, 14, 28) & 
				ModelRuns$Efficacy_CamelControls == unique(ModelRuns$Efficacy_CamelControls)[1], 	PrintCols]

#### Maximum impact of reactive (no camels)
# cases (mean, lower and upper)
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.0, 	"PropCasesAverted_Mean"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.0, 	"PropCasesAverted_LowerCrI"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.0, 	"PropCasesAverted_UpperCrI"])

# deaths (mean, lower and upper)
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.0, 	"DeathsAverted_CF_Mean"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.0, 	"DeathsAverted_CF_LowerCrI"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.0, 	"DeathsAverted_CF_UpperCrI"])


#### Maximum impact of reactive (with 30% camels)
# cases (mean, lower and upper)
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.3, 	"PropCasesAverted_Mean"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.3, 	"PropCasesAverted_LowerCrI"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.3, 	"PropCasesAverted_UpperCrI"])

# deaths (mean, lower and upper)
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.3, 	"DeathsAverted_CF_Mean"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.3, 	"DeathsAverted_CF_LowerCrI"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.3, 	"DeathsAverted_CF_UpperCrI"])


#### Maximum impact of reactive (with 50% camels)
# cases (mean, lower and upper)
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.5, 	"PropCasesAverted_Mean"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.5, 	"PropCasesAverted_LowerCrI"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.5, 	"PropCasesAverted_UpperCrI"])

# deaths (mean, lower and upper)
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.5, 	"DeathsAverted_CF_Mean"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.5, 	"DeathsAverted_CF_LowerCrI"])
max(ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$Efficacy_CamelControls == 0.5, 	"DeathsAverted_CF_UpperCrI"])



ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & 
				ModelRuns$ReactLevel == "REGIONAL" &
				ModelRuns$Efficacy_Start == 1 & 
				ModelRuns$ImplementationDelay %in% c(0, 14) & 
				ModelRuns$Efficacy_CamelControls == unique(ModelRuns$Efficacy_CamelControls)[1], 	PrintCols]




