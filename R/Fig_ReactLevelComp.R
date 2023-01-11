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
ModelRuns = read.table(file = file.path(CppOutputDirectory, paste0("ModelRunsSummary.txt")), header = T, sep = "\t") 

Efficacy_CamelControls = 0
ModelRuns_Reactive_Hosp = ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$ReactLevel == "HOSPITAL" 	 & 
				ModelRuns$ImmunityDelay == 14 & ModelRuns$Efficacy_CamelControls == Efficacy_CamelControls,]
ModelRuns_Reactive_Region = ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$ReactLevel  == "REGIONAL"	 & 
				ModelRuns$ImmunityDelay == 14 & ModelRuns$Efficacy_CamelControls == Efficacy_CamelControls,]
ModelRuns_Reactive_Nation = ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$ReactLevel  == "NATIONAL"	 & 
				ModelRuns$ImmunityDelay == 14 & ModelRuns$Efficacy_CamelControls == Efficacy_CamelControls,]

CF_quantity = "PropCasesAverted_Mean"

## On average, what is the difference between hospital and regional reactive policies, irrespective of efficacy and implementation delay
mean(ModelRuns_Reactive_Region[, CF_quantity] / ModelRuns_Reactive_Hosp[, CF_quantity])
mean(ModelRuns_Reactive_Nation[, CF_quantity] / ModelRuns_Reactive_Hosp[, CF_quantity])

RatioOrDifference = function(a, b, Ratio_bool = TRUE) if (Ratio_bool) a / b else a - b

Delay_Cols 		= bpy.colors(n = length(unique(ModelRuns$ImplementationDelay)))
LWD				= 3
XaxisVec 		= seq(5, 100, 5)
Xlab 			= "Efficacy (%)"

Ratio_bool 	= TRUE
CF_quantity = "PropCasesAverted_Mean"

for (CF_quantity in c("PropCasesAverted_Mean", "DeathsAverted_CF_Mean"))
{
	if (CF_quantity == "PropCasesAverted_Mean") CF_Name = "Cases averted" else CF_Name = "Deaths averted"  
	png(filename = file.path(here("Plots"), paste0("National_vs_Regional_vs_Hospital_Reactive_", CF_Name, ".png")), res = 500, units = "in", width = 12, height = 7.5)
	LayoutMatrix = matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE)
	LayoutMatrix = cbind(LayoutMatrix, c(7,7))
	layout(LayoutMatrix, widths = c(1,1,1,0.4))
	
	for (Ratio_bool in c(TRUE, FALSE))
	{
		if (Ratio_bool) Ratio_Char = "/" 		else Ratio_Char = "-"
		if (Ratio_bool) Ratio_Name = "ratio" 	else Ratio_Name = "difference"
		
		# Get y-axis ranges so that each comparison plotted on same scale. Ranges irrespective of efficacy and implementation delay
		Ranges = c()
		for (WhichComparison in 1:3)
		{
			if (WhichComparison == 1) { DFrame_1 = ModelRuns_Reactive_Hosp		; 	DFrame_2 = ModelRuns_Reactive_Region } ## Regional vs. Hospital
			if (WhichComparison == 2) { DFrame_1 = ModelRuns_Reactive_Hosp		;	DFrame_2 = ModelRuns_Reactive_Nation } ## National vs. Hospital
			if (WhichComparison == 3) { DFrame_1 = ModelRuns_Reactive_Region	; 	DFrame_2 = ModelRuns_Reactive_Nation } ## National vs. Regional
			
			# overall range (i.e. irrespective of efficacy and implementation delay)
			Ranges = c(Ranges, range(RatioOrDifference(DFrame_2[, CF_quantity], DFrame_1[, CF_quantity], Ratio_bool)))
		}
		MinYaxis = min(Ranges)
		MaxYaxis = max(Ranges)
		
		for (WhichComparison in 1:3)
		{
			if (WhichComparison == 1) { DFrame_1 = ModelRuns_Reactive_Hosp		; 	DFrame_2 = ModelRuns_Reactive_Region } ## Regional vs. Hospital
			if (WhichComparison == 2) { DFrame_1 = ModelRuns_Reactive_Hosp		;	DFrame_2 = ModelRuns_Reactive_Nation } ## National vs. Hospital
			if (WhichComparison == 3) { DFrame_1 = ModelRuns_Reactive_Region	; 	DFrame_2 = ModelRuns_Reactive_Nation } ## National vs. Regional
			
			OrigMAR = par("mar")
			if (WhichComparison == 1) Main = paste0("Regional vs. Hospital")
			if (WhichComparison == 2) Main = paste0("National vs. Hospital")
			if (WhichComparison == 3) Main = paste0("National vs. Regional")
			
			if (Ratio_bool) Main = paste0(Main, "\nRatio of proportions ("		, CF_Name, ")")  	else 
				Main = paste0(Main, "\nDifference in proportions ("	, CF_Name, ")")
			
			DelayIndex = 1
			for (DelayIndex in 1:length(unique(ModelRuns$ImplementationDelay)))
			{
				Delay 		= unique(ModelRuns$ImplementationDelay)[DelayIndex]
				YaxisVec 	= RatioOrDifference(	
						DFrame_2	[which(DFrame_2$ImplementationDelay == Delay), CF_quantity], 
						DFrame_1	[which(DFrame_1$ImplementationDelay == Delay), CF_quantity], Ratio_bool)
				
				if (DelayIndex == 1)
					plot	(XaxisVec, YaxisVec, col = Delay_Cols[DelayIndex], ylim = c(MinYaxis, MaxYaxis), 
							type = "b", lwd = LWD, cex.main = 1.5,
							main = Main, xlab = Xlab, ylab = Ratio_Name, cex.lab = 1.6) else 
					lines	(XaxisVec, YaxisVec, col = Delay_Cols[DelayIndex], lwd = LWD, type = "b")
			}
		}
	}
	par(mar = rep(0, 4))
	plot(NULL , xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
	legend("center", legend = paste0(unique(ModelRuns$ImplementationDelay), " days"), col = Delay_Cols, lwd = LWD,
			pch = NA, lty = 1, cex = 1.6, bty = "n")
	par(mar = OrigMAR) ### reset default margins
	dev.off()
}




