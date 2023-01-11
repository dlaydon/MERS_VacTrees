# Libraries
library(dplyr)
library(tidyr)
library(viridis)
require(sp)

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

Hill = function(Time, HillPower, Halflife)
{
	ToverH_RaisedToPower = (Time/Halflife)^HillPower
	return(ToverH_RaisedToPower / (1 + ToverH_RaisedToPower));
}

## Make waning picture
PlotPolygon = function(XValues, LowerLine, UpperLine, MeanLine = NULL, PolyCOL = "pink", MiddleCol = "red", ALPHA = 0.5, LWD = 4, type = "l")
{
	polygon	(c(XValues,rev(XValues)),c(LowerLine,rev(UpperLine)), col = PolyCOL, border = NA)
	if (!is.null(MeanLine)) lines(x = XValues, y = MeanLine, col = MiddleCol, lwd = LWD, type = "b")
}

HillPower = 4
ParDefaults = par()
OrigMAR 	= ParDefaults$mar

############################################################################
OrigMAI = par ("mai") ### good to record if using layout functions
OrigMAR = par ("mar") ### good to record if using layout functions

ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
PlotsDirectory		= file.path(ProjectDirectory, "Plots"		)

source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))
source(file.path(R_ScriptDirectory, "LayoutPlottingFunctions.R"))

# Import summary model runs
ModelRuns = read.table(file = file.path(ProjectDirectory, paste0("ModelRunsSummary.txt")), header = T, sep = "\t")
ModelRuns = ModelRuns[	ModelRuns$MinDay_Pruning == 0 	& 
				ModelRuns$MaxDay_Pruning == 548 , ] 

DurationNames = paste0(unique(ModelRuns$VaccineDuration), " year")
DurationNames[which(DurationNames != "1 year")] = paste0(DurationNames[which(DurationNames != "1 year")], "s")

png(filename = file.path(here("Plots"), paste0("WaningPic2.png")), res = 200, units = "in", width = 13, height = 5.5)
layout(matrix(1:3, nrow = 1, ncol = 3), widths = c(1,1,0.3))
par(xpd = TRUE)
# Make model runs copy without efficacy = 0 or infinite duration.
MRCopy = ModelRuns[which(ModelRuns$Efficacy_Start != 0 & ModelRuns$VaccineDuration != Inf),]

LWD = 4
CF_quantity = "PropCasesAverted_Mean"
XVac = seq(0,20,0.05)

CEXMAIN = 1.7
CEXLAB = 1.5
CEXAXIS  = 1.4

## (A) Exponential waning
plot(NA, xlim = c(0, 20), ylim = c(0,100), ylab = "Percentage efficacy remaining", xlab = "Years post dose", 
		main = "Exponential waning\nby mean duration of protection", cex.lab = CEXLAB, cex.main = CEXMAIN, cex.axis = CEXAXIS)
DurIndex = 2
Cols 		= bpy.colors(n = length(unique(MRCopy$VaccineDuration)))
AlphaCols	= bpy.colors(n = length(unique(MRCopy$VaccineDuration)), alpha = 1)
Cols == unique(Cols)
for (DurIndex in 1:length(unique(MRCopy$VaccineDuration))) 
	PlotPolygon(XVac, rep(0, length(XVac)), UpperLine = 100 * exp(-XVac/unique(MRCopy$VaccineDuration)[DurIndex]), PolyCOL = AlphaCols[DurIndex])
text(x = -2.5, y = 112, label = "A", cex = 3.5)


## (B) Hill (sigmoidal) waning

plot(NA, xlim = c(0, 20), ylim = c(0,100), ylab = "Percentage efficacy remaining", xlab = "Years post dose", 
		main = "Sigmoidal waning\nby half-life", cex.lab = CEXLAB, cex.main = CEXMAIN, cex.axis = CEXAXIS)
DurIndex = 2
Cols 		= bpy.colors(n = length(unique(MRCopy$VaccineDuration)))
AlphaCols	= bpy.colors(n = length(unique(MRCopy$VaccineDuration)), alpha = 1)
Cols == unique(Cols)
for (DurIndex in 1:length(unique(MRCopy$VaccineDuration))) 
	PlotPolygon(XVac, rep(0, length(XVac)), UpperLine = 100 * (1 - Hill(XVac, HillPower, Halflife = unique(MRCopy$VaccineDuration)[DurIndex])), PolyCOL = AlphaCols[DurIndex])
text(x = -2.5, y = 112, label = "B", cex = 3.5)

par(mar = rep(0, 4), xpd = FALSE)
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
LegCols = bpy.colors(n = length(unique(MRCopy$VaccineDuration)))
legend("center", title = "Mean duration/\nHalf-life", legend = DurationNames[-1], 
		col = LegCols, lwd = 8, pch = NA, lty = 1, cex = 1.5, bty = "n")
par(mar = OrigMAR) ### reset default margins
dev.off()




## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## Make camel only cases/deaths averted picture
ModelRunsCamelsOnly = ModelRuns[ModelRuns$Efficacy_Start == 0, ]

PlotPolygon = function(XValues, LowerLine, UpperLine, MeanLine = NULL, PolyCOL = "pink", MiddleCol = "red", ALPHA = 0.5, LWD = 4, type = "l")
{
	polygon	(c(XValues,rev(XValues)),c(LowerLine,rev(UpperLine)), col = PolyCOL, border = NA)
	if (!is.null(MeanLine)) lines(x = XValues, y = MeanLine, col = MiddleCol, lwd = LWD, type = "b")
}

png(file.path(PlotsDirectory, "CamelsOnlyCasesDeathsAverted.png"), units = "in", width = 11, height = 5, res = 500)
par(mfrow = c(1,2))
par(xpd = TRUE)

CEXAXIS = 1.2
CEXLAB 	= 1.3
plot(NA, xlim = c(0.1, max(ModelRunsCamelsOnly$Efficacy_CamelControls)), 
		ylim = 0:1, #xaxt = "n", yaxt = "n", 
		bty = "n", 
		main = "Cases averted",
		xlab = "Camel control measure effectiveness", 
		ylab = "Proportion cases averted", cex.axis = CEXAXIS, cex.lab = CEXLAB)
text(x = 0.03, y = 1.23, label = "A", cex = 4)
PlotPolygon(XValues = ModelRunsCamelsOnly$Efficacy_CamelControls, 
		LowerLine = ModelRunsCamelsOnly$PropCasesAverted_LowerCrI, 
		UpperLine = ModelRunsCamelsOnly$PropCasesAverted_UpperCrI,
		MeanLine  = ModelRunsCamelsOnly$PropCasesAverted_Mean	)
plot(NA, xlim = c(0.1, max(ModelRunsCamelsOnly$Efficacy_CamelControls)), 
		ylim = 0:1, #xaxt = "n", yaxt = "n", 
		bty = "n", 
		main = "Deaths averted",
		xlab = "Camel control measure effectiveness", 
		ylab = "Proportion deaths averted", cex.axis = CEXAXIS, cex.lab = CEXLAB)
text(x = 0.03, y = 1.23, label = "B", cex = 4)
PlotPolygon(XValues = ModelRunsCamelsOnly$Efficacy_CamelControls, 
		LowerLine = ModelRunsCamelsOnly$DeathsAverted_CF_LowerCrI, 
		UpperLine = ModelRunsCamelsOnly$DeathsAverted_CF_UpperCrI,
		MeanLine  = ModelRunsCamelsOnly$DeathsAverted_CF_Mean	)
dev.off()







