
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
OrigMAI = par ("mai") ### record if using layout functions
OrigMAR = par ("mar") ### record if using layout functions

ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
PlotsDirectory		= file.path(ProjectDirectory, "Plots"		)

source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))
source(file.path(R_ScriptDirectory, "LayoutPlottingFunctions.R"))

# Import summary model runs
ModelRuns = read.table(file = file.path(CppOutputDirectory, paste0("ModelRunsSummary.txt")), header = T, sep = "\t") 

ListOfCounterfactuals = c("PropCasesAverted", "PropCasesAverted_HCW", "PropCasesAverted_nHCW", 
		"DeathsAverted_CF", "DeathsAverted_CF_HCW", "DeathsAverted_CF_nHCW", "FinalCaseDate_CF", "Peak", "Peak_Reduction")
ZVar_Chars 	= c(ListOfCounterfactuals, "Peak_Reduction")


#### Residual efficacy
PlotNum = 1
Stat = "_Mean"
VariableInQuestion 	= paste0(ZVar_Chars[PlotNum], Stat)

ModelRunsCopy = ModelRuns
ModelRunsCopy$Efficacy_Start = ModelRunsCopy$Efficacy_Start * 100

ProactiveIndices = which(ModelRunsCopy$VacCampStrategy == "PROACTIVE")


ModelRunsCopy$ResidEff 						= ModelRunsCopy$Efficacy_Start
ModelRunsCopy$ResidEff[ProactiveIndices] 	= 
		(ModelRunsCopy$Efficacy_Start * exp(-ModelRunsCopy$TimeSinceVaccination / ModelRunsCopy$VaccineDuration))[ProactiveIndices]

ModelRunsCopy$ResidEff = 
		(ModelRunsCopy$Efficacy_Start * exp(-ModelRunsCopy$TimeSinceVaccination / ModelRunsCopy$VaccineDuration))
ModelRunsCopy$ResidEff[-ProactiveIndices] = ModelRunsCopy$Efficacy_Start[-ProactiveIndices]


Efficacy_CamelControls = 0.0
# Proactive indices
P_Ind = which(ModelRunsCopy$VacCampStrategy == "PROACTIVE" & ModelRunsCopy$Efficacy_CamelControls == Efficacy_CamelControls)

ReactLevel 	= "HOSPITAL"
#ReactLevel 	= "REGIONAL"

png(filename = here("Plots/ResidEffsCasesAverted.png"), units = "in", res = 500, width = 11, height = 4)

LayoutMatrix = matrix(1:3, nrow = 1)
LayoutMatrix = cbind(LayoutMatrix, max(LayoutMatrix) + 1) # legend
LayoutMatrix = rbind(LayoutMatrix, max(LayoutMatrix) + 1) # shared x-axis
LayoutMatrix = cbind(max(LayoutMatrix) + 1, LayoutMatrix) # shared y-axis
layout(LayoutMatrix, widths = c(0.06, 1,1,1,0.7), heights = c(1, 0.05))

COLS = c("black", "grey", "red", "green", "darkblue", "orange", "purple")
Cols 		= COLS[-c(1,2)]

for (ReactLevel in unique(ModelRunsCopy$ReactLevel))
{
	plot(ModelRunsCopy$ResidEff[P_Ind], ModelRunsCopy[P_Ind, VariableInQuestion],
			ylim = c(0, max(ModelRunsCopy[ModelRunsCopy$Efficacy_CamelControls == Efficacy_CamelControls, VariableInQuestion])),
			main = ReactLevel,
			#ylim = 0:1, 
			col = COLS[1], lwd = 4, cex.main = 1.7,
			xlab = "", ylab = "")
	
	points(ModelRunsCopy$Efficacy_Start[P_Ind], ModelRunsCopy[P_Ind, VariableInQuestion], 
			col = COLS[2], lwd = 1)
	
	arrows(x0 = ModelRunsCopy$Efficacy_Start[P_Ind], y0 = ModelRunsCopy[P_Ind, VariableInQuestion], 
			x1 = ModelRunsCopy$ResidEff[P_Ind], y1 = ModelRunsCopy[P_Ind, VariableInQuestion], 
			length = 0, col = "grey")
	
	ImpDelays 	= c(0,8,14,22,28)
	for (ImpDelayIndex in 1:length(ImpDelays))
	{
		ImpDelay = ImpDelays[ImpDelayIndex]
		
		# Reactive indices for this react time/ implementation delay.
		R_Ind 	= which(ModelRunsCopy$VacCampStrategy == "REACTIVE" & ModelRunsCopy$Efficacy_CamelControls == Efficacy_CamelControls & 
						ModelRunsCopy$ReactLevel == ReactLevel & ModelRunsCopy$ImplementationDelay == ImpDelay)
		
		points(ModelRunsCopy$Efficacy_Start[R_Ind], ModelRunsCopy[R_Ind, VariableInQuestion], 
				col = Cols[ImpDelayIndex], lwd = 3, 
#				pch = ImpDelayIndex + 1, 
				pch = 1, 
				cex = 3)
	}
}
CEXAXIS = 1.8

par(mar = rep(0,4))
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend("left", legend = c(paste0(ImpDelays, "-day react time"), "Proactive (residual)", "Proactive (initial)"), 
		col = c(Cols, "black", "grey"), cex = 1.5, lty = NA, 
#		pch = 1,
		pch = c(1:length(ImpDelays) + 1, 1, 1),
		lwd = c(rep(3, length(ImpDelays)),4,1), pt.cex)
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
text(0.41, 0.55, "Efficacy (%)", cex = CEXAXIS)
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
text(0.5, 0.54, "Proportion of cases averted", srt = 90, cex = CEXAXIS)
par(mar = OrigMAR) ### reset default margins

dev.off()

