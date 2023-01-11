# Libraries
library(dplyr)
library(tidyr)
library(viridis)

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
require(stringr)

ParDefaults = par()
OrigMAR 	= ParDefaults$mar


############################################################################
OrigMAI = par ("mai") ### record if using layout functions
OrigMAR = par ("mar") ### record if using layout functions

ProjectDirectory 	= here()
RawDataDirectory 	= file.path(ProjectDirectory, "Data"		)
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
PlotsDirectory		= file.path(ProjectDirectory, "Plots"		)

source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))
source(file.path(R_ScriptDirectory, "LayoutPlottingFunctions.R"))

# Import processed data for Cpp
DATA 		= read.table(file = file.path(RawDataDirectory, "MERS_forCpp.txt"), header = TRUE)
Day_0  		= as.Date		("2013-01-01") #### 1st Jan 2013 is Day zero. 
Day_Final  	= Day_0 + max(DATA$onset) # note Different from Day_Final in MakeLineListForCpp.R script.

TotalDeaths 		= sum(DATA$Dead)
TotalDeaths_HCW 	= length(which(DATA$HCW == 1 & DATA$Dead == 1))
TotalDeaths_nHCW 	= length(which(DATA$HCW == 0 & DATA$Dead == 1))
LastDay_Data_Int 	= max(DATA$onset)

## make epidemic curve for Data
DATA_EpiCurves 			= rep(NA, LastDay_Data_Int + 1)
DATA_EpiCurves_Deaths 	= rep(NA, LastDay_Data_Int + 1)
for (Day in 0:LastDay_Data_Int) DATA_EpiCurves			[Day + 1] = length(which(DATA$onset == Day))
for (Day in 0:LastDay_Data_Int) DATA_EpiCurves_Deaths	[Day + 1] = length(which(DATA$onset == Day & DATA$Dead == 1))

ConvertDailyToWeeklyInc = function(DailyInc)
{
	DaysToConsider 	= seq(1, length(DailyInc), by = 7)
	WeeklyInc 		= c(0, diff(cumsum(DailyInc), lag = 7))[DaysToConsider] # have to include zero for the first week 
	if (sum(WeeklyInc, na.rm = T) != sum(DailyInc, na.rm = T)) 
		WeeklyInc[length(WeeklyInc)] = sum(DailyInc, na.rm = T) - sum(WeeklyInc, na.rm = T)
	return(WeeklyInc)
}
DATA_EpiCurves_Weekly 			= ConvertDailyToWeeklyInc(DATA_EpiCurves)
DATA_EpiCurves_Weekly_Deaths	= ConvertDailyToWeeklyInc(DATA_EpiCurves_Deaths)

cbind(DATA_EpiCurves_Weekly, DATA_EpiCurves_Weekly_Deaths)
DATA_EpidemicPeakSize_Weekly 	= max(DATA_EpiCurves_Weekly)

DateLabels 			= c("Jan '13", "Apr '13", "Jul '13", "Oct '13", "Jan '14", "Apr '14", "Jul '14")
DatesForAxis 		= as.Date(c("2013-01-01", "2013-04-01", "2013-07-01", "2013-10-01", "2014-01-01", "2014-04-01", "2014-07-01"))
DatesForAxis_Int 	= round(as.numeric(DatesForAxis - Day_0) / 7)
DatesForAxis_Int[1] = 1 # not zero 

PlotPolygon = function(XValues, LowerLine, UpperLine, MeanLine = NULL, PolyCOL = "pink", MiddleCol = "red", ALPHA = 0.5, LWD = 4, type = "l")
{
	polygon	(c(XValues,rev(XValues)),c(LowerLine,rev(UpperLine)), col= PolyCOL, border = NA)
	if (!is.null(MeanLine)) lines(x = XValues, y = MeanLine, col = MiddleCol, lwd = LWD, type = "l")
}
CasesCol = "grey70"              
DeathsCol = "grey70"             
DeathsCol = "black"             
HCW_Col = "blue"               
nHCW_Col = "green2"            
Deaths_HCW_Col = "gold"      
Deaths_nHCW_Col = "orangered"      
CEXLEG = 1.17
LWD = 2


################ =========================== ################ =========================== ################ =========================== 
################ =========================== National epi curves with counterfactuals superimposed

# Import summary model runs
ModelRuns = read.table(file = file.path(CppOutputDirectory, paste0("ModelRunsSummary.txt")), header = T, sep = "\t") 
ModelRuns = ModelRuns[ModelRuns$MinDay_Pruning== 0 & ModelRuns$MaxDay_Pruning == 548,]

## Decide on subset to plot

# Proactive, ExpWaning, NoWindow, 70% Efficacy, 10yr Duration, 8 year wait, No Camels
ModelRunsSubsetToPlot_1 = ModelRuns[	ModelRuns$VacCampStrategy 				== "REACTIVE" 					& 
				ModelRuns$WaningInReactive				== 1							&
				ModelRuns$ReactLevel 					%in% c("HOSPITAL", "NATIONAL")	&
				ModelRuns$ImmunityDelay 				== 14 							& 
				ModelRuns$VaccineDuration 				== 10 							& 
				ModelRuns$ImplementationDelay 			%in% c(0, 28)					&
				ModelRuns$Efficacy_Start 				%in% c(0.5,0.95) 				& 
				ModelRuns$ExpWaning 					== 1 							& 
				ModelRuns$Efficacy_CamelControls 		== 0.0							,] 
ModelRunsSubsetToPlot_1 = ModelRunsSubsetToPlot_1[c(1,2,5,6,3,4,7,8),]

ModelRunsSubsetToPlot_2 = ModelRuns[	ModelRuns$VacCampStrategy 				== "PROACTIVE" 		& 
				ModelRuns$ImmunityDelay 				== 0 				& 
				ModelRuns$VaccineDuration 				%in% c(2,20)		& 
				ModelRuns$TimeSinceVaccination 			%in% c(1, 10) 		& 
				ModelRuns$Efficacy_Start 				%in% c(0.5,0.95)	& 
				ModelRuns$ExpWaning 					== 1 				& 
				ModelRuns$Efficacy_CamelControls 		== 0.0				,]
ModelRunsSubsetToPlot_2 = ModelRunsSubsetToPlot_2[order(ModelRunsSubsetToPlot_2$Efficacy_Start),]
ModelRunsSubsetToPlot = rbind(ModelRunsSubsetToPlot_1, ModelRunsSubsetToPlot_2)					

par (xpd = TRUE, mfcol = c(2,4))

png(filename = file.path(here("Plots/"), paste0("ExampleEpiCurves.png")), res = 100, units = "in", width = 10, height = 10)

PlotCounter = 2 ## Meta rownames
PlotCounter = PlotCounter + 4 # rownames
FirstPlotMatrix 	= matrix(1:8 + PlotCounter, nrow = 2)
SecondPlotMatrix 	= matrix(1:8 + max(FirstPlotMatrix), nrow = 2)

PlotMatrix = rbind(FirstPlotMatrix, SecondPlotMatrix)
PlotMatrix = cbind(c(1,1,2,2), 3:6, PlotMatrix)
PlotMatrix = PlotMatrix + 2 ## Meta Colnames
PlotMatrix = rbind(rep(NA, 6), PlotMatrix)
PlotMatrix[1,3:6] = c(1,1,2,2) ## Meta Colnames
PlotMatrix[1,1:2] = max(PlotMatrix, na.rm = T) + 1

layout(PlotMatrix, widths = c(0.7,0.5,2,2,2,2), heights = c(1,3,3,3,3))
par(mar = c(0,0,0,0))
# MetaColnames
MetaColnames = c("50% efficacy", "95% efficacy")
for (MetaColname in 1:2)
{
	plot(NA, xlim = 0:1, ylim = 0:1, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
	text(0.5, 0.5, MetaColnames[MetaColname], cex = 3, font = 2) ### don't know why you have to reverse the rownames, but it otherwise the names come out in the wrong order. 
}
# Meta Rownames
MetaRownames = c("Reactive", "Proactive")
for (MetaRowname in 1:2)
{
	plot(NA, xlim = 0:1, ylim = 0:1, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
	text(0.5, 0.5, MetaRownames[MetaRowname], cex = 3, font = 2, srt = 90) ### don't know why you have to reverse the rownames, but it otherwise the names come out in the wrong order. 
}
#Rownames
Rownames = c("Hospital-level", "National-level", "20-year mean duration", "2-year mean duration")
for (Row in 1:4)
{
	plot(NA, xlim = 0:1, ylim = 0:1, xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
	text(0.5, 0.5, Rownames[Row], cex = 1.6, font = 2, srt = 90) ### don't know why you have to reverse the rownames, but it otherwise the names come out in the wrong order. 
}
#par(mar = OrigMAR, xpd = TRUE)
par(mar = c(4.1, 2.1, 3.1, 2.1), xpd = TRUE)
for (MR_Index in 1:dim(ModelRunsSubsetToPlot)[1])
{
	MR = ModelRunsSubsetToPlot[MR_Index,]
	
	OutputString 			= paste0("_", MR$OutputFolderNames)
	EpiCurvesFileName 		= file.path(CppOutputDirectory, paste0("CF_EpiCurves" 			, OutputString, ".txt"))
	CF_EpiCurves 			= read.table(file = EpiCurvesFileName		, header = T, sep = "\t")
	
	if (MR$VacCampStrategy == "REACTIVE") 
		Title = paste0(MR$ImplementationDelay	, "-day react time") else
		Title = paste0(MR$TimeSinceVaccination	, "-year wait")
	
	b <- barplot(DATA_EpiCurves_Weekly, col = CasesCol, border = NA, width = 1,	main = Title, cex.main = 1.6)
	PlotPolygon(XValues = b[,1], LowerLine = CF_EpiCurves$LowerCrI, 
			UpperLine = CF_EpiCurves$UpperCrI, PolyCOL = "pink", MiddleCol = "red", 
			MeanLine = CF_EpiCurves$Mean, LWD = LWD)
	text(x = b[DatesForAxis_Int,1], y = -25, label = DateLabels, srt = 35) # have to place text at coords where barplot is actually drawn - because for some reason bar widths aren't equal to 1.
	axis(1, at=b[DatesForAxis_Int,1], labels=FALSE)
	legend("topleft", legend = c("Data", "Mean", "CrI"), 
			col = c("grey70", "red", "pink"), lty = 1, lwd = 4)
}
dev.off()


################ =========================== ################ =========================== ################ =========================== 
################ =========================== Hospital and regional epi curves, with reaction windows superimposed.


## rearrange and reorder Data so that first region has first case etc.
DATA 		= read.table(file = file.path(RawDataDirectory, "MERS_forCpp.txt"), header = TRUE)

DATA = DATA[order(DATA$onset), ]
RegionalData 	= read.csv	(file = file.path(RawDataDirectory, "RegionalHealthcareData_DJL.txt")	, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#DATA$Region 	= match(DATA$Region, unique(DATA$Region))
#DATA$Hospital 	= match(DATA$Hospital, unique(DATA$Hospital))

data.frame(Orig = DATA$Hospital, Reordered = match(DATA$Hospital, unique(DATA$Hospital)))

ImmunityDelay = 14

HospitalOrRegionChar = "Hospital"
HospitalOrRegionChar = "Region"
for (HospitalOrRegionChar in c("Hospital", "Region"))
{
	if (HospitalOrRegionChar == "Region")
		png(filename = file.path(here("Plots/"), paste0("EpiCurvesBy", HospitalOrRegionChar, ".png")), units = "in", width = 12, height = 10, res = 300)
	if (HospitalOrRegionChar == "Hospital")
		png(filename = file.path(here("Plots/"), paste0("EpiCurvesBy", HospitalOrRegionChar, ".png")), units = "in", width = 24, height = 15, res = 300)
	if (HospitalOrRegionChar == "Region")	par(mfrow = c(3,4), xpd = FALSE) else
	if (HospitalOrRegionChar == "Hospital")	par(mfrow = c(6,8), xpd = FALSE, mar = c(2.1, 4.1, 4.1, 2.1))
	
	Onset_Week 		= floor(DATA$onset / 7)
	AllCases_Inc 	= rep(NA, max(Onset_Week)); 
	Cases_HCW 		= rep(NA, max(Onset_Week)); 
	
	HospitalOrRegion = 3
	for (HospitalOrRegion in 1:length(unique(DATA[, HospitalOrRegionChar])))
	{
		if (HospitalOrRegionChar == "Region") 
		{
			RegionName = unique(RegionalData$RegionNameInData[which(RegionalData$RegionNumber == (Region - 1))])
			RegionName = RegionName[1]
		}
		
		# clear incidence vectors
		AllCases_Inc[1:length(AllCases_Inc)	] = NA
		Cases_HCW	[1:length(Cases_HCW)	] = NA
		
		DATAdummy 		= DATA[DATA[, HospitalOrRegionChar] == (HospitalOrRegion - 1),]
		FirstCaseTime 	= min(DATAdummy$onset)
		HCW_CaseTimes 	= DATAdummy$onset[DATAdummy$HCW == 1]
		Onset_Week 		= floor(DATAdummy$onset / 7)

		AllCases_Inc	[as.numeric(names(table(Onset_Week)))] 								= table(Onset_Week)
		Cases_HCW		[as.numeric(names(table(Onset_Week[which(DATAdummy$HCW == 1)])))] 	= table(Onset_Week[which(DATAdummy$HCW == 1)])
		
		if (HospitalOrRegionChar == "Hospital" & length(which(DATAdummy$HCW == 1)) == 0) next 
		
		if (HospitalOrRegionChar == "Region") 	
			YCoordForDateLabels = max(AllCases_Inc, na.rm = TRUE) * -0.2 else 
			YCoordForDateLabels = max(AllCases_Inc, na.rm = TRUE) * -0.1 
		
		
		PlotTitle = paste0(HospitalOrRegionChar, " ", HospitalOrRegion)
		if (HospitalOrRegionChar == "Region") PlotTitle = paste0(PlotTitle, ": ", RegionName)
		b <- barplot(AllCases_Inc, col = CasesCol, border = NA, width = 1, main = PlotTitle)
		
		## plot vertical lines
		abline(v = 1.2 * (seq(0,28,2) + FirstCaseTime + ImmunityDelay)/7, col = "lightgreen", alpha = 0.5) # 1 = bar width, 0.2 = gap. Hence 1.2
		
		par(xpd = TRUE)
		text(x = b[DatesForAxis_Int,1], y = YCoordForDateLabels, label = DateLabels, srt = 35) # have to place text at coords where barplot is actually drawn - because for some reason bar widths aren't equal to 1.
		axis(1, at = b[DatesForAxis_Int,1], labels = FALSE)
		if (HospitalOrRegionChar == "Region") 
			legend("topleft", legend = c(paste0("All cases (n = ", dim(DATAdummy)[1], " cases)"), paste0("HCW cases (n = ", length(which(!is.na(Cases_HCW))), " cases)")), col = c(CasesCol, HCW_Col), pch = 15, cex = CEXLEG)
		if (length(which(DATAdummy$HCW == 1)) > 0) 
			barplot(Cases_HCW, col = HCW_Col, border = NA, add = TRUE, width = 1)
		par(xpd = FALSE)
	}
	dev.off()
}















