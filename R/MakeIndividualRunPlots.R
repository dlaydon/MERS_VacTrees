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

# import model runs
ModelRuns = read.table(file = file.path(CppOutputDirectory, paste0("ModelRunsSummary.txt")), header = T, sep = "\t") 

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

DidRunFail 	= rep(0, dim(ModelRuns)[1])

for (MR_index in 1:dim(ModelRuns)[1])
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
	
	
	ModelRunPlot_dir = file.path(PlotsDirectory, OutputSubDir) ### given only one file, perhaps better to have all the files together, not subdivided by folder. 
	dir.create(ModelRunPlot_dir, showWarnings = FALSE)
	
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

	
	#### clean up. 
	FilenamesCompleted = c(FilenamesCompleted, OutputString)
	rm(CHAINS, CF_CHAINS, OutputSubDir, NumRows, NumCols, NumChainPlotsPerPNG, CF_EpiCurves, CF_EpiCurves_Deaths)
	CloseOpenPlotDevices()	
	cat("\n")
}

