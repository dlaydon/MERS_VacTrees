
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

ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
PlotsDirectory		= file.path(ProjectDirectory, "Plots"		)

source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))
source(file.path(R_ScriptDirectory, "LayoutPlottingFunctions.R"))

COLNAMES_TransContribMat_Props 	= c("Contrib_SameHosp", "Contrib_SameRegion", "Contrib_DiffRegion", "Contrib_Reservoir")
COLNAMES_TransContribMat_Abs 	= c("Contrib_SameHosp_MeanAbs", "Contrib_SameRegion_MeanAbs", "Contrib_DiffRegion_MeanAbs", "Contrib_Reservoir_MeanAbs")
COLNAMES_TransContribMat 		= c(COLNAMES_TransContribMat_Props, COLNAMES_TransContribMat_Abs)

# Import summary model runs
ModelRuns = read.table(file = file.path(CppOutputDirectory, paste0("ModelRunsSummary.txt")), header = T, sep = "\t") 
TransmissionContribMatrix = ModelRuns[, COLNAMES_TransContribMat]

BarCols = c("orange", "green", "blue", "red") ## choose non-horrible colours later.
LWD = 5
CEXAXIS = 1
CEXLAB = 1.5
CEXMAIN = 1.4

TransmissionBarplot = function(Indices, PlotTitle, Xaxis = NULL, Xlab = "Efficacy (%)", 
		COLNAMES = COLNAMES_TransContribMat_Props, YLIM = 0.61)
{
	TransMatDummy = TransmissionContribMatrix[Indices, COLNAMES]
	if (is.null(Xaxis)) Xaxis = ModelRuns$Efficacy_Start[Indices] * 100
	
	plot(NA, type = "l", col = BarCols[1], xlim = range(Xaxis), ylim = c(0, YLIM), lwd = LWD, main = PlotTitle, 
			ylab = "", xlab = Xlab, cex.axis = CEXAXIS, cex.lab = CEXLAB, cex.main = CEXMAIN)
	for (Type in c(1,2,4,3))
		lines(Xaxis, t(TransMatDummy)[Type,], col = BarCols[Type], lwd = LWD)
}

ImpDelay = 28

IndicesHOSPTIAL = 	
		ModelRuns$VacCampStrategy 			== "REACTIVE" 			& 
		ModelRuns$ReactLevel 				== "HOSPITAL" 			& 
		ModelRuns$Efficacy_CamelControls 	== 0					&			
		ModelRuns$ImmunityDelay 			== 14					&			
		ModelRuns$ImplementationDelay 		== ImpDelay					
IndicesREGIONAL = 	ModelRuns$VacCampStrategy 			== "REACTIVE" 			& 
		ModelRuns$ReactLevel 				== "REGIONAL" 			& 
		ModelRuns$Efficacy_CamelControls 	== 0					&			
		ModelRuns$ImmunityDelay 			== 14					&			
		ModelRuns$ImplementationDelay 		== ImpDelay					
IndicesNATIONAL = 	ModelRuns$VacCampStrategy 	== "REACTIVE" 			& 
		ModelRuns$ReactLevel 				== "NATIONAL" 		& 
		ModelRuns$Efficacy_CamelControls 	== 0				&			
		ModelRuns$ImmunityDelay 			== 14				&			
		ModelRuns$ImplementationDelay 		== 2				
IndicesNATIONAL_wCamels = 	ModelRuns$VacCampStrategy 	== "REACTIVE" 			& 
		ModelRuns$ReactLevel 				== "NATIONAL" 		& 
		ModelRuns$Efficacy_CamelControls 	== 0.3				&			
		ModelRuns$ImmunityDelay 			== 14				&			
		ModelRuns$ImplementationDelay 		== ImpDelay				

IndicesPROACTIVE_noCamels = ModelRuns$VacCampStrategy 	== "PROACTIVE" 	& 
		ModelRuns$Efficacy_CamelControls 	== 0					&	
		ModelRuns$VaccineDuration 	== 5				&			
		ModelRuns$TimeSinceVaccination 			== 0.5	
IndicesPROACTIVE_wCamels = ModelRuns$VacCampStrategy 	== "PROACTIVE" 	& 
		ModelRuns$Efficacy_CamelControls 	== 0.3					&	
		ModelRuns$VaccineDuration 	== 5				&			
		ModelRuns$TimeSinceVaccination 			== 0.5	
IndicesCAMELS = ModelRuns$Efficacy_Start 	== 0	 		


png(filename = here("Plots/TransmissionContribs.png"), units = "in", res = 500, width = 10, height = 6.5)
LayoutMatrix = matrix(1:8, nrow = 2, ncol = 4, byrow = TRUE)
LayoutMatrix = cbind(max(LayoutMatrix) + 1, LayoutMatrix)
layout(LayoutMatrix, widths = c(0.1,1,1,1,1))

TransmissionBarplot(IndicesHOSPTIAL			, "(A) Reactive\n(Hospital)")
TransmissionBarplot(IndicesREGIONAL			, "(B) Reactive\n(Regional)")
TransmissionBarplot(IndicesNATIONAL			, "(C) Reactive\n(National)")
TransmissionBarplot(IndicesNATIONAL_wCamels	, "(D) Reactive (National)\n30% effective camel controls")
TransmissionBarplot(IndicesCAMELS, "(E) Camels controls only", Xaxis = ModelRuns$Efficacy_CamelControls[IndicesCAMELS] * 100, Xlab = "Effectiveness (%)")
TransmissionBarplot(IndicesPROACTIVE_noCamels, "(F) Proactive")
TransmissionBarplot(IndicesPROACTIVE_wCamels, "(G) Proactive (with 30%\neffective camel controls)")

par(mar = rep(0, 4))
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend("center", legend = c("Within hospital", "Within region", "Between regions", "Animal reservoir"), 
		title = "Transmission type",
		col = BarCols, pch = 15, cex = 2, bty = "n", pt.cex = 4)
plot(NA, xlim = 0:1, ylim = c(0,1), yaxs='i', xaxs='i', xaxt = "n", yaxt = "n", bty="n")
text(0.5, 0.5, "Transmission contribution (proportional)", srt = 90, cex = 2)
par(mar = OrigMAR) ### reset default margins
dev.off()

png(filename = here("Plots/TransmissionContribsAbs.png"), units = "in", res = 500, width = 10, height = 6.5)
LayoutMatrix = matrix(1:8, nrow = 2, ncol = 4, byrow = TRUE)
LayoutMatrix = cbind(max(LayoutMatrix) + 1, LayoutMatrix)
layout(LayoutMatrix, widths = c(0.1,1,1,1,1))

TransmissionBarplot(IndicesHOSPTIAL				, "(A) Reactive\n(Hospital)"							, COLNAMES = COLNAMES_TransContribMat_Abs, YLIM = 380)
TransmissionBarplot(IndicesREGIONAL				, "(B) Reactive\n(Regional)"							, COLNAMES = COLNAMES_TransContribMat_Abs, YLIM = 380)
TransmissionBarplot(IndicesNATIONAL				, "(C) Reactive\n(National)"							, COLNAMES = COLNAMES_TransContribMat_Abs, YLIM = 380)
TransmissionBarplot(IndicesNATIONAL_wCamels		, "(D) Reactive (National)\n30% effective camel controls", COLNAMES = COLNAMES_TransContribMat_Abs, YLIM = 380)
TransmissionBarplot(IndicesCAMELS				, "(E) Camels controls only"							, COLNAMES = COLNAMES_TransContribMat_Abs, YLIM = 380,
		Xaxis = ModelRuns$Efficacy_CamelControls[IndicesCAMELS] * 100, Xlab = "Effectiveness (%)")
TransmissionBarplot(IndicesPROACTIVE_noCamels	, "(F) Proactive"										, COLNAMES = COLNAMES_TransContribMat_Abs, YLIM = 380)
TransmissionBarplot(IndicesPROACTIVE_wCamels	, "(G) Proactive (with 30%\neffective camel controls)"	, COLNAMES = COLNAMES_TransContribMat_Abs, YLIM = 380)

par(mar = rep(0, 4))
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend("center", legend = c("Within hospital", "Within region", "Between regions", "Animal reservoir"), 
		title = "Transmission type",
		col = BarCols, pch = 15, cex = 2, bty = "n", pt.cex = 4)
plot(NA, xlim = 0:1, ylim = c(0,1), yaxs='i', xaxs='i', xaxt = "n", yaxt = "n", bty="n")
text(0.5, 0.5, "Transmission contribution (absolute)", srt = 90, cex = 2)
par(mar = OrigMAR) ### reset default margins
dev.off()


