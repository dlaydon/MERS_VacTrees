

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
require(sp)


ParDefaults = par()
OrigMAR 	= ParDefaults$mar

############################################################################
#options(width = 172L)
options(width = 108L)
OrigMAI = par ("mai") ### good to record if using layout functions
OrigMAR = par ("mar") ### good to record if using layout functions


ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
PlotsDirectory		= file.path(ProjectDirectory, "Plots"		)

source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))
source(file.path(R_ScriptDirectory, "LayoutPlottingFunctions.R"))


# Import summary model runs
ModelRuns = read.table(file = file.path(CppOutputDirectory, paste0("ModelRunsSummary.txt")), header = T, sep = "\t") 
DurationNames = paste0(unique(ModelRuns$VaccineDuration), " year")
DurationNames[which(DurationNames != "1 year")] = paste0(DurationNames[which(DurationNames != "1 year")], "s")



## Make waning picture
PlotPolygon = function(XValues, LowerLine, UpperLine, MeanLine = NULL, PolyCOL = "pink", MiddleCol = "red", ALPHA = 0.5, LWD = 4, type = "l")
{
	polygon	(c(XValues,rev(XValues)),c(LowerLine,rev(UpperLine)), col = PolyCOL, border = NA)
	if (!is.null(MeanLine)) lines(x = XValues, y = MeanLine, col = MiddleCol, lwd = LWD, type = "b")
}
png(filename = file.path(here("Plots"), paste0("WaningPic2.png")), res = 200, units = "in", width = 7, height = 5.5)
layout(matrix(1:2, nrow = 1, ncol = 2), widths = c(1,0.3))
# Make model runs copy without efficacy = 0 or infinite duration.
MRCopy = ModelRuns[which(ModelRuns$Efficacy_Start != 0 & ModelRuns$VaccineDuration != Inf),]

LWD = 4
CF_quantity = "PropCasesAverted_Mean"
XVac = seq(0,20,0.05)
plot(NA, xlim = c(0, 20), ylim = c(0,100), ylab = "Percentage efficacy remaining", xlab = "Years post dose", 
		main = "Percentage initial vaccine efficacy\nremaining over time by mean duration", cex.lab = 1.4, cex.main = 1.4, cex.axis = 1.2)
DurIndex = 2
Cols 		= bpy.colors(n = length(unique(MRCopy$VaccineDuration)))
AlphaCols	= bpy.colors(n = length(unique(MRCopy$VaccineDuration)), alpha = 1)
Cols == unique(Cols)
for (DurIndex in 1:length(unique(MRCopy$VaccineDuration))) 
	PlotPolygon(XVac, rep(0, length(XVac)), UpperLine = 100 * exp(-XVac/unique(MRCopy$VaccineDuration)[DurIndex]), PolyCOL = AlphaCols[DurIndex])

par(mar = rep(0, 4))
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
LegCols = bpy.colors(n = length(unique(MRCopy$VaccineDuration)))
legend("left", title = "Mean vaccine\nduration", legend = DurationNames[-1], col = LegCols, lwd = 8, pch = NA, lty = 1, cex = 1.1, bty = "n")
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
CEXAXIS = 1.2
CEXLAB 	= 1.3
plot(NA, xlim = c(0.1, max(ModelRunsCamelsOnly$Efficacy_CamelControls)), 
		ylim = 0:1, #xaxt = "n", yaxt = "n", 
		bty = "n", 
		main = "(A) Proportion cases averted by\ncamel control measure effectiveness",
		xlab = "Camel control measure effectiveness", 
		ylab = "Proportion cases averted", cex.axis = CEXAXIS, cex.lab = CEXLAB)
PlotPolygon(XValues = ModelRunsCamelsOnly$Efficacy_CamelControls, 
		LowerLine = ModelRunsCamelsOnly$PropCasesAverted_LowerCrI, 
		UpperLine = ModelRunsCamelsOnly$PropCasesAverted_UpperCrI,
		MeanLine  = ModelRunsCamelsOnly$PropCasesAverted_Mean	)
plot(NA, xlim = c(0.1, max(ModelRunsCamelsOnly$Efficacy_CamelControls)), 
		ylim = 0:1, #xaxt = "n", yaxt = "n", 
		bty = "n", 
		main = "(B) Proportion deaths averted by\ncamel control measure effectiveness",
		xlab = "Camel control measure effectiveness", 
		ylab = "Proportion deaths averted", cex.axis = CEXAXIS, cex.lab = CEXLAB)
PlotPolygon(XValues = ModelRunsCamelsOnly$Efficacy_CamelControls, 
		LowerLine = ModelRunsCamelsOnly$DeathsAverted_CF_LowerCrI, 
		UpperLine = ModelRunsCamelsOnly$DeathsAverted_CF_UpperCrI,
		MeanLine  = ModelRunsCamelsOnly$DeathsAverted_CF_Mean	)
dev.off()







