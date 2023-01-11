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
require(sp)
require(stringr)

ParDefaults = par()
OrigMAR 	= ParDefaults$mar

############################################################################
OrigMAI = par ("mai") ### record if using layout functions
OrigMAR = par ("mar") ### record if using layout functions

ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
PlotsDirectory		= file.path(ProjectDirectory, "Plots"		)

source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))
source(file.path(R_ScriptDirectory, "LayoutPlottingFunctions.R"))

# Import summary model runs
ModelRuns = read.table(file = file.path(ProjectDirectory, paste0("ModelRunsSummary.txt")), header = T, sep = "\t") 
str(ModelRuns)

PlotNames = c("Direct", "Indirect")
ListOfCounterfactuals = c("Cases_DirectlyProtected"	, "Cases_IndirectlyProtected")
LegendNames	= c("Proportion of\ncases averted", "Proportion of\ncases averted")
ZVar_Chars 	= ListOfCounterfactuals

# subset to only relevant model runs
ModelRuns = ModelRuns[ModelRuns$Efficacy_CamelControls == 0, ]
ModelRuns = ModelRuns[ModelRuns$ExpWaning == 1, ]
ModelRuns = ModelRuns[ModelRuns$MinDay_Pruning == 0, ]
ModelRuns = ModelRuns[ModelRuns$MaxDay_Pruning == 548, ]
ModelRuns = ModelRuns[!is.na(ModelRuns$Cases_DirectlyProtected_Mean), ]

ModelRuns$Proportion = ModelRuns$Cases_DirectlyProtected_Mean / ModelRuns$PropCasesAverted_Mean

all(ModelRuns$PropCasesAverted_Mean == ModelRuns$Cases_DirectlyProtected_Mean + ModelRuns$Cases_IndirectlyProtected_Mean)
mean(abs(ModelRuns$PropCasesAverted_Mean - (ModelRuns$Cases_DirectlyProtected_Mean + ModelRuns$Cases_IndirectlyProtected_Mean))/ModelRuns$PropCasesAverted_Mean)

#par(mfcol = c(4,2))
png(filename = file.path(here("Plots/DirectContrib.png")), units = "in", res = 300, width = 12, height = 8)
layout(cbind(matrix(1:8, nrow = 2, byrow = TRUE), c(9,9)))
Duration = 20

CEXLAB = 1.3
CEXAXIS = 1.3
EffSub = c(0.2,0.4, 0.6, 0.8, 1)
Cols = bpy.colors(n = length(unique(EffSub)) + 1)[-1]
LineTypes = c(1, 2,3)

for (Duration in c(1, 5, 10, 20))
{
	if (Duration == 1)	YearOrYears = "year" else YearOrYears = "years"
	
	plot(NA, xlim = c(0, 28), ylim = c(0,0.4), xlab = "React time (days)", ylab = "Proportion", cex.lab = CEXLAB, cex.axis = CEXAXIS,
			main = paste0("Reactive\nDuration of protection = ", Duration, " ", YearOrYears), font.main = 1)
	TypeCounter = 1

	for (ReactLevel in c("HOSPITAL", "REGIONAL", "NATIONAL"))
	{
		ColCounter = 1
		for (Eff in EffSub)
		{
			ModelRuns_Reactive = ModelRuns[ModelRuns$VacCampStrategy 	== "REACTIVE" 				& 
							ModelRuns$ReactLevel 					== ReactLevel				& 
							ModelRuns$WaningInReactive				== 1			&
							ModelRuns$VaccineDuration				== Duration			&
							ModelRuns$Efficacy_Start 				== Eff						& 
							ModelRuns$ImmunityDelay 				== 14						, ]
			dim(ModelRuns_Reactive)
			ModelRuns_Reactive[, 1:15]
			
			ModelRuns_Reactive$Cases_DirectlyProtected_Mean
			Proportions = ModelRuns_Reactive$Cases_DirectlyProtected_Mean / ModelRuns_Reactive$PropCasesAverted_Mean
			lines(unique(ModelRuns_Reactive$ImplementationDelay), Proportions, col = Cols[ColCounter], lwd = 5, lty = LineTypes[TypeCounter])
			
			ColCounter = ColCounter + 1
		}
		TypeCounter = TypeCounter + 1
		
	}
}

for (Duration in c(1, 5, 10, 20))
{
	if (Duration == 1)	YearOrYears = "year" else YearOrYears = "years"
	
	plot(NA, xlim = c(0, 10), ylim = c(0,0.4), xlab = "Wait (years)", ylab = "Proportion", cex.lab = CEXLAB, cex.axis = CEXAXIS,
			main = paste0("Proactive\nDuration of protection = ", Duration, " ", YearOrYears), font.main = 1)
	TypeCounter = 1
	
	ColCounter = 1
	for (Eff in EffSub)
	{
		ModelRuns_Proactive = ModelRuns[ModelRuns$VacCampStrategy 	== "PROACTIVE" 	& 
						ModelRuns$VaccineDuration				== Duration			&
						ModelRuns$Efficacy_Start 				== Eff				& 
						ModelRuns$ImmunityDelay 				== 0				, ]
		dim(ModelRuns_Proactive)
		ModelRuns_Proactive[, 1:15]
		
		ModelRuns_Proactive$Cases_DirectlyProtected_Mean
		Proportions = ModelRuns_Proactive$Cases_DirectlyProtected_Mean / ModelRuns_Proactive$PropCasesAverted_Mean
		lines(unique(ModelRuns_Proactive$TimeSinceVaccination), Proportions, col = Cols[ColCounter], lwd = 5, lty = LineTypes[TypeCounter])
		
		ColCounter = ColCounter + 1
	}
	TypeCounter = TypeCounter + 1
}
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend("top", title = "React level", legend = c("Hospital", "Regional", "National"), lty = LineTypes, lwd = 5, bty = "n", cex = 2)
legend("bottom", , title = "Efficacy", legend = paste0(EffSub * 100, "%"), lty = 1, col = Cols[1:length(EffSub)], lwd = 5, bty = "n", cex = 2)
par(mar = OrigMAR) ### reset default margins
dev.off()








