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
OrigMAI = par ("mai") ### good to record if using layout functions
OrigMAR = par ("mar") ### good to record if using layout functions

ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppRootDirectory 	= file.path(ProjectDirectory, "MERS_Vac"	, "MERS_Vac") 
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
RawDataDirectory 	= file.path(ProjectDirectory, "Data"		)
source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))

# Import processed data for Cpp
DATA 		= read.table(file = file.path(RawDataDirectory, "MERS_forCpp.txt"), header = TRUE)
DATA_safe 	= DATA

Onset_Week 		= floor(DATA$onset / 7)
NumWeeks		= max(Onset_Week)
AllCases_Inc 	= rep(NA, max(Onset_Week)); 
AllCases_Inc[as.numeric(names(table(Onset_Week)))] = table(Onset_Week)

str(DATA)
attach(DATA)
Day_0  				= as.Date		("2013-01-01") #### 1st Jan 2013 is Day zero. 
Day_Final  			= Day_0 + max(DATA$onset) # note Different from Day_Final in MakeLineListForCpp.R script.
LastDay_Data_Int 	= max(DATA$onset)

FindWaitTime  = function(NumberCasesToHit, TimePeriod, Collated = Collated_National)
{
	Diff = Collated[, TimePeriod]
	WaitTime = min(which(Diff >= NumberCasesToHit))
	if (WaitTime == Inf) WaitTime = NA
	return(WaitTime)
}
FindWaitTimes = function(NumberCasesToHit, Collated = Collated_National)
{
	TimePeriods = c("Day", "Week", "TwoWeeks", "Month", "TwoMonths")
	WaitTimes = rep(NA, length(TimePeriods))
	names(WaitTimes) = TimePeriods
	for (TimePeriod in TimePeriods)
		WaitTimes[TimePeriod] = FindWaitTime(NumberCasesToHit, TimePeriod = TimePeriod, Collated = Collated)
	WaitTimes = c(NumberCasesToHit, WaitTimes)
	names(WaitTimes)[1] = "Trigger"
	return(WaitTimes)
}

### === NATIONAL VERSION OF FIGURE
# Import processed data for Cpp
DATA 		= read.table(file = file.path(RawDataDirectory, "MERS_forCpp.txt"), header = TRUE)
DATA_safe 	= DATA

str(DATA)
attach(DATA)
Day_0  				= as.Date		("2013-01-01") #### 1st Jan 2013 is Day zero. 
Day_Final  			= Day_0 + max(DATA$onset) # note Different from Day_Final in MakeLineListForCpp.R script.
LastDay_Data_Int 	= max(DATA$onset)

DATA_EpiCurves 		= rep(NA, LastDay_Data_Int + 1)
for (Day in 0:LastDay_Data_Int) DATA_EpiCurves[Day + 1] = length(which(DATA$onset == Day))

WeeklyDiff 		= c(cumsum(DATA_EpiCurves)[1:7 ], diff(cumsum(DATA_EpiCurves), lag = 7 ))
BiWeeklyDiff 	= c(cumsum(DATA_EpiCurves)[1:14], diff(cumsum(DATA_EpiCurves), lag = 14))
MonthlyDiff 	= c(cumsum(DATA_EpiCurves)[1:30], diff(cumsum(DATA_EpiCurves), lag = 30))
BiMonthlyDiff 	= c(cumsum(DATA_EpiCurves)[1:60], diff(cumsum(DATA_EpiCurves), lag = 60))

Collated_National = data.frame(Date = Day_0 + (1:length(DATA_EpiCurves) - 1), 
		Day = DATA_EpiCurves, Week = WeeklyDiff, TwoWeeks = BiWeeklyDiff, Month = MonthlyDiff, TwoMonths = BiMonthlyDiff)

FindWaitTime = function(NumberCasesToHit, TimePeriod, Collated = Collated_National)
{
	Diff = Collated[, TimePeriod]
	WaitTime = min(which(Diff >= NumberCasesToHit))
	if (WaitTime == Inf) WaitTime = NA
	return(WaitTime)
}
FindWaitTimes = function(NumberCasesToHit, Collated = Collated_National)
{
	TimePeriods = c("Day", "Week", "TwoWeeks", "Month", "TwoMonths")
	WaitTimes = rep(NA, length(TimePeriods))
	names(WaitTimes) = TimePeriods
	for (TimePeriod in TimePeriods)
		WaitTimes[TimePeriod] = FindWaitTime(NumberCasesToHit, TimePeriod = TimePeriod, Collated = Collated)
	WaitTimes = c(NumberCasesToHit, WaitTimes)
	names(WaitTimes)[1] = "Trigger"
	return(WaitTimes)
}

TriggerDFrame = FindWaitTimes(NumberCasesToHit = 1)
for (Trigger in 2:30) TriggerDFrame = rbind(TriggerDFrame, FindWaitTimes(NumberCasesToHit = Trigger))
rownames(TriggerDFrame) = NULL
TriggerDFrame = as.data.frame(TriggerDFrame)
TriggerDFrame


Onset_Week 		= floor(DATA$onset / 7)
AllCases_Inc 	= rep(NA, max(Onset_Week)); 
AllCases_Inc[as.numeric(names(table(Onset_Week)))] = table(Onset_Week)

DateLabels 			= c("Jan '13", "Apr '13", "Jul '13", "Oct '13", "Jan '14", "Apr '14", "Jul '14")
DatesForAxis 		= as.Date(c("2013-01-01", "2013-04-01", "2013-07-01", "2013-10-01", "2014-01-01", "2014-04-01", "2014-07-01"))
DatesForAxis_Int 	= round(as.numeric(DatesForAxis - Day_0) / 7)
DatesForAxis_Int[1] = 1 # not zero 
DayOrDays = c("day", rep("days", 4))

png(filename = file.path(here("Plots/TriggersNational_AndEpiCurve.png")), units = "in", res = 300, width = 8, height = 5)
layout(matrix(1:2, nrow = 1), widths = c(1, 0.3))
par(xpd = TRUE, mar = c(5.1,4.1,4.1,5.1))
b <- barplot(AllCases_Inc, col = "grey", border = NA, width = 1, main = paste0(""), ylab = "Cases")
text(x = b[DatesForAxis_Int,1], y = -10, label = DateLabels, srt = 35) # have to place text at coords where barplot is actually drawn - because for some reason bar widths aren't equal to 1.
axis(1, at = b[DatesForAxis_Int, 1], labels = FALSE)
par(new = TRUE)
PlotTitle = "National"
plot(NA, xlim = c(0, length(AllCases_Inc)), ylim = c(0, max(TriggerDFrame$Trigger)), type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side = 4, at = pretty(range(TriggerDFrame$Trigger)))
TriggerDFrame[is.na(TriggerDFrame)] = LastDay_Data_Int + 10
TrigCols = adjustcolor(palette()[1:dim(TriggerDFrame)[2]], alpha.f = 0.7)
for (Trigger in TriggerDFrame$Trigger)
	for (Col in 2:dim(TriggerDFrame)[2])
	{
		segments(0, Trigger, TriggerDFrame[Trigger, Col]/7, Trigger, col = TrigCols[Col], lwd = 2)
		if (TriggerDFrame[Trigger, Col] != LastDay_Data_Int + 10)
			points(TriggerDFrame[Trigger, Col]/7, Trigger, pch = 4, col = TrigCols[Col], lwd = 1.2)
	}
mtext("Trigger value", side = 4, line = 3)
warnings()
par(mar = rep(0, 4), xpd = FALSE)
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend("center", title = "Timeframe", legend = paste(c(1,7,14,30,60), DayOrDays),  
		col = TrigCols[2:dim(TriggerDFrame)[2]], lwd = 2, pch = 4, lty = 1, cex = 1.2, bty = "n")
par(mar = OrigMAR) ### reset default margins
dev.off()


### === REGIONAL VERSION OF FIGURE


DATA 			= read.table(file = file.path(RawDataDirectory, "MERS_forCpp.txt"), header = TRUE)
RegionalData 	= read.csv	(file = file.path(RawDataDirectory, "RegionalHealthcareData_DJL.txt")	, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

png(filename = file.path(here("Plots/TriggersRegional_AndEpiCurve.png")), units = "in", res = 500, width = 12, height = 10)
LayoutMatrix = matrix(1:12, nrow = 4, ncol = 3, byrow = T)
LayoutMatrix
layout(LayoutMatrix)

DATA 			= DATA[order(DATA$onset), ]
data.frame(1:length(unique(DATA$Region)), unique(DATA$Region))

head(DATA, 20)

Region = 1
for (Region in 1:length(unique(DATA_safe$Region)))
{
	SubData = DATA[DATA$Region == (Region - 1),]
	
	RegionName = unique(RegionalData$RegionNameInData[which(RegionalData$RegionNumber == (Region - 1))])
	
	DATA_EpiCurves 		= rep(NA, LastDay_Data_Int + 1)
	for (Day in 0:LastDay_Data_Int) DATA_EpiCurves[Day + 1] = length(which(SubData$onset == Day))
	
	WeeklyDiff 		= c(cumsum(DATA_EpiCurves)[1:7 ], diff(cumsum(DATA_EpiCurves), lag = 7 ))
	BiWeeklyDiff 	= c(cumsum(DATA_EpiCurves)[1:14], diff(cumsum(DATA_EpiCurves), lag = 14))
	MonthlyDiff 	= c(cumsum(DATA_EpiCurves)[1:30], diff(cumsum(DATA_EpiCurves), lag = 30))
	BiMonthlyDiff 	= c(cumsum(DATA_EpiCurves)[1:60], diff(cumsum(DATA_EpiCurves), lag = 60))
	
	Collated_National = data.frame(Date = Day_0 + (1:length(DATA_EpiCurves) - 1), 
			Day = DATA_EpiCurves, Week = WeeklyDiff, TwoWeeks = BiWeeklyDiff, Month = MonthlyDiff, TwoMonths = BiMonthlyDiff)
	
	TriggerDFrame = FindWaitTimes(NumberCasesToHit = 1)
	for (Trigger in 2:30) TriggerDFrame = rbind(TriggerDFrame, FindWaitTimes(NumberCasesToHit = Trigger))
	rownames(TriggerDFrame) = NULL
	TriggerDFrame = as.data.frame(TriggerDFrame)
	
	
	### instead plot full epidemic curve - with vertical lines denoting date that a particular trigger is reached.
	### have multiple panels for multiple triggers or multiple time periods.
	
	Onset_Week 		= floor(SubData$onset / 7)
	AllCases_Inc 	= rep(NA, NumWeeks); 
	AllCases_Inc[as.numeric(names(table(Onset_Week)))] = table(Onset_Week)
	
	DateLabels 			= c("Jan '13", "Apr '13", "Jul '13", "Oct '13", "Jan '14", "Apr '14", "Jul '14")
	DatesForAxis 		= as.Date(c("2013-01-01", "2013-04-01", "2013-07-01", "2013-10-01", "2014-01-01", "2014-04-01", "2014-07-01"))
	DatesForAxis_Int 	= round(as.numeric(DatesForAxis - Day_0) / 7)
	DatesForAxis_Int[1] = 1 # not zero 
	DayOrDays = c("day", rep("days", 4))
	
	par(xpd = TRUE, mar = c(5.1,4.1,4.1,5.1))
	YCoordForDateLabels = max(AllCases_Inc, na.rm = TRUE) * -0.2
	b <- barplot(AllCases_Inc, col = "grey", border = NA, width = 1, main = paste0(""), ylab = "Cases", cex.lab = 1.3)
	text(x = b[DatesForAxis_Int,1], y = YCoordForDateLabels, label = DateLabels, srt = 35, cex = 1.12) # have to place text at coords where barplot is actually drawn - because for some reason bar widths aren't equal to 1.
	axis(1, at = b[DatesForAxis_Int, 1], labels = FALSE)
	par(new = TRUE)
	PlotTitle = paste0("Region ", Region, ": ", RegionName)
	plot(NA, xlim = c(0, length(AllCases_Inc)), ylim = c(0, max(TriggerDFrame$Trigger)), 
			main = PlotTitle,  font.main = 1, cex.main = 2, 
			type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "")
	axis(side = 4, at = pretty(range(TriggerDFrame$Trigger)))
	TriggerDFrame[is.na(TriggerDFrame)] = LastDay_Data_Int + 10
	TrigCols = adjustcolor(palette()[1:dim(TriggerDFrame)[2]], alpha.f = 0.3)
	for (Trigger in TriggerDFrame$Trigger)
		for (Col in 2:dim(TriggerDFrame)[2])
		{
			segments(0, Trigger, TriggerDFrame[Trigger, Col]/7, Trigger, col = TrigCols[Col], lwd = 1)
			if (TriggerDFrame[Trigger, Col] != LastDay_Data_Int + 10)
				points(TriggerDFrame[Trigger, Col]/7, Trigger, pch = 4, col = TrigCols[Col], lwd = 1.2)
		}
	mtext("Trigger", side = 4, line = 2.5, cex = 0.8)
	warnings()
}
par(mar = rep(0, 4), xpd = FALSE)
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend("center", title = "Timeframe", legend = paste(c(1,7,14,30,60), DayOrDays),  
		col = TrigCols[2:dim(TriggerDFrame)[2]], lwd = 3, pch = 4, lty = 1, cex = 2.5, bty = "n")
par(mar = OrigMAR) ### reset default margins
dev.off()




DATA = read.table(file = file.path(RawDataDirectory, "MERS_forCpp.txt"), header = TRUE)

## Next figure should show, for given trigger and time period, 
## how many regions and hospitals would even catch epidemic.

Trigger = 4
Region = 5 

## rearrange and reorder Data so that first region has first case etc.
DATA 			= DATA[order(DATA$onset), ]

MetaTrigDFrame 	= expand.grid(
		Wait = c(30, 180, 365, 550), 
		ReactLevel = c("Hospital", "Region"), 
		Timeframe = c("Day", "Week", "TwoWeeks", "Month", "TwoMonths"), 
		Trigger = 1:20)

MetaTrigDFrame$Proportion = rep(0, dim(MetaTrigDFrame)[1]) # initialize all to zero, then increment in loop below
head(MetaTrigDFrame, 20)

HospitalOrRegionChar = "Hospital"
HospitalOrRegionChar = "Region"
for (HospitalOrRegionChar in c("Hospital", "Region"))
{
	Onset_Week 		= floor(DATA$onset / 7)
	AllCases_Inc 	= rep(NA, max(Onset_Week)); 
	Cases_HCW 		= rep(NA, max(Onset_Week)); 
	
	HospitalOrRegion = 1
	for (HospitalOrRegion in 1:length(unique(DATA[, HospitalOrRegionChar])))
	{
		DATAdummy = DATA[DATA[, HospitalOrRegionChar] == (HospitalOrRegion),]
		
		DATA_EpiCurves_RegionOrHospital = rep(NA, LastDay_Data_Int + 1)
		for (Day in 0:LastDay_Data_Int) DATA_EpiCurves_RegionOrHospital[Day + 1] = length(which(DATAdummy$onset == Day))
		
		WeeklyDiff 		= c(cumsum(DATA_EpiCurves_RegionOrHospital)[1:7 ]	, diff(cumsum(DATA_EpiCurves_RegionOrHospital), lag = 7 ))
		BiWeeklyDiff 	= c(cumsum(DATA_EpiCurves_RegionOrHospital)[1:14]	, diff(cumsum(DATA_EpiCurves_RegionOrHospital), lag = 14))
		MonthlyDiff 	= c(cumsum(DATA_EpiCurves_RegionOrHospital)[1:30]	, diff(cumsum(DATA_EpiCurves_RegionOrHospital), lag = 30))
		BiMonthlyDiff 	= c(cumsum(DATA_EpiCurves_RegionOrHospital)[1:60]	, diff(cumsum(DATA_EpiCurves_RegionOrHospital), lag = 60))
		
		Collated_RegionOrHospital  = data.frame(Date = Day_0 + (1:length(DATA_EpiCurves_RegionOrHospital) - 1),
				Day = DATA_EpiCurves_RegionOrHospital, Week = WeeklyDiff, TwoWeeks = BiWeeklyDiff, Month = MonthlyDiff, TwoMonths = BiMonthlyDiff)
		
		TriggerDFrame = FindWaitTimes(NumberCasesToHit = 1, Collated = Collated_RegionOrHospital)
		for (Trigger in 2:30) 
			TriggerDFrame = rbind(TriggerDFrame, FindWaitTimes(NumberCasesToHit = Trigger, Collated = Collated_RegionOrHospital))
		rownames(TriggerDFrame) = NULL
		TriggerDFrame = as.data.frame(TriggerDFrame)
		TriggerDFrame[is.na(TriggerDFrame)] = Inf
		
		for (Wait in unique(MetaTrigDFrame$Wait))
			for (Trigger in unique(MetaTrigDFrame$Trigger))
				for (Timeframe in unique(MetaTrigDFrame$Timeframe))
				{
					# Find index within MetaTrigDFrame
					Index = which(	MetaTrigDFrame$Wait 		== Wait 				& 
									MetaTrigDFrame$Timeframe 	== Timeframe 			& 
									MetaTrigDFrame$ReactLevel 	== HospitalOrRegionChar & 
									MetaTrigDFrame$Trigger 		== Trigger 				)

					# if trigger (for this trigger threshold and Timeframe) was reached before Wait, increment.
					if (TriggerDFrame[Trigger, Timeframe] <= Wait + min(DATA$onset))
						MetaTrigDFrame[Index, "Proportion"] = MetaTrigDFrame[Index, "Proportion"] + 1
				}
	}
}
MetaTrigDFrame[MetaTrigDFrame$ReactLevel == "Hospital"	, "Proportion"] = MetaTrigDFrame[MetaTrigDFrame$ReactLevel == "Hospital", "Proportion"] / length(unique(DATA$Hospital	))
MetaTrigDFrame[MetaTrigDFrame$ReactLevel == "Region"	, "Proportion"] = MetaTrigDFrame[MetaTrigDFrame$ReactLevel == "Region"	, "Proportion"] / length(unique(DATA$Region		))


for (Wait in unique(MetaTrigDFrame$Wait))
	for (ReactLevel in unique(MetaTrigDFrame$ReactLevel))
	{
		MetaTrigDFrame_Sub = MetaTrigDFrame[MetaTrigDFrame$ReactLevel == ReactLevel & MetaTrigDFrame$Wait == Wait, ]
		
		p <- ggplot(MetaTrigDFrame_Sub, aes(x = Trigger, y = Proportion)) + 
				geom_line(aes(colour = Timeframe), size  = 1) +
				geom_point(aes(colour = Timeframe, shape = Timeframe), size  = 1.5) +
				labs(y = paste0(ReactLevel, "s"), x = "") + 
				ylim(0,1) + 
				theme_minimal() + 
#				theme(legend.position = "none") + 
				ggtitle(paste0(ReactLevel, "-level, within ", Wait, " days")) 
		GPlotWLegend <- p + guides(fill = guide_legend(title = "Timeframe")) 
		legend <- get_legend(GPlotWLegend)
		assign(paste0("Legend_Wait_", Wait),  legend)	
		
		p  <- p + theme(legend.position = "none")
		assign(paste0("Wait_", Wait, "_", ReactLevel, "_Plot"), p)
	}

Filename = paste0("TriggersReachedByRegionAndHospital")
png(filename = file.path(here("Plots/"), paste0(Filename, ".png")),	res = 500, units = "in", width = 9, height = 11)
#pdf(file = file.path(here("Plots/"), paste0(Filename, ".pdf")),	width = 21, height = 17)

LayoutMatrix = matrix(1:(length(unique(MetaTrigDFrame$Wait)) * length(unique(MetaTrigDFrame$ReactLevel))), 
		nrow = length(unique(MetaTrigDFrame$Wait)), byrow = T)
LayoutMatrix = cbind(LayoutMatrix, max(LayoutMatrix) + 1)

grid.arrange(
		
		Wait_30_Hospital_Plot	, 
		Wait_30_Region_Plot		, 
		Wait_180_Hospital_Plot	, 
		Wait_180_Region_Plot	, 
		Wait_365_Hospital_Plot	, 
		Wait_365_Region_Plot	, 
		Wait_550_Hospital_Plot	, 
		Wait_550_Region_Plot	, 
		
		Legend_Wait_30			,
		
		layout_matrix = LayoutMatrix,
		widths 	= c(1,1,0.5)			,
		
		nrow 	= length(unique(MetaTrigDFrame$Wait)), 
		bottom 	= textGrob("Trigger threshold (MERS-CoV infections)"		, gp = gpar(fontsize = 15, font = 1), hjust = 0.7), 
		left 	= textGrob("Proportion hospitals/regions where trigger reached"			, gp = gpar(fontsize = 15, font = 1), rot = 90)) ## would be better to use grid arrange with characters but leave for now.
dev.off()







