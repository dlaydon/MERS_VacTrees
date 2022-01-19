
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

XVar_Char = "Efficacy_Start"
YVar_Char = "VaccineDuration"

ZVar_Chars 	= c(ListOfCounterfactuals, "Peak_Reduction")
LegendNames	= c("Proportion of\ncases averted"	, "Proportion of\ncases averted\n(healthcare\nworkers only)" , "Proportion of\ncases averted\n(excluding\nhealthcare workers)", 
		"Proportion of\ndeaths averted"	, "Proportion of\ndeaths averted\n(healthcare\nworkers only)", "Proportion of\ndeaths averted\n(excluding\nhealthcare workers)", 
		"Final\ncase day"				, "Epidemic\npeak size", "Proportional\nreduction in\nepidemic peak")
str(ModelRuns$VaccineDuration)
str(ModelRuns)

# Make copy of model runs so can alter variables. E.g. non-uniform spacing between values in TimeSinceVaccination, or 0 being coded as no waning (i.e. infinite duration).

ModelRunsCopy = ModelRuns
ModelRunsCopy = ModelRunsCopy[ModelRunsCopy$TimeSinceVaccination 	%in% c(0.5,1,2,4,8), ]
ModelRunsCopy = ModelRunsCopy[ModelRunsCopy$Efficacy_CamelControls 	%in% unique(ModelRunsCopy$Efficacy_CamelControls)[c(1,4,6)], ]
### need ModelRunsCopy$VaccineDuration as factor
ModelRunsCopy$VaccineDuration[which(ModelRunsCopy$VaccineDuration == 0)] = Inf
ModelRunsCopy$Efficacy_Start 			= ModelRunsCopy$Efficacy_Start * 100
ModelRunsCopy$Efficacy_CamelControls 	= ModelRunsCopy$Efficacy_CamelControls * 100
ModelRunsCopy$VaccineDuration 			= as.factor(ModelRunsCopy$VaccineDuration		)
ModelRunsCopy$TimeSinceVaccination 		= as.factor(ModelRunsCopy$TimeSinceVaccination	)
TimesSinceVaccinationNames 				= paste0(unique(ModelRunsCopy$TimeSinceVaccination), " year")
Efficacy 								= unique(ModelRunsCopy$Efficacy_Start)[18]

#ModelRuns_subset[ModelRuns_subset$VaccineDuration == 10, ]
#any(is.na(ModelRuns_subset))

PlotNum = 1
for (PlotNum in 1:length(ZVar_Chars))
{
	cat(paste0("PlotNum ", PlotNum, ": ", ZVar_Chars[PlotNum], "\n"))
	
	Stat = "_Mean"
	for (Stat in c("_Mean", "_LowerCrI", "_UpperCrI"))
	{
		cat(paste0(Stat, "\n"))
		
		EffCamelIndex = 1
		for (EffCamelIndex in 1:length(unique(ModelRunsCopy$Efficacy_CamelControls)))
		{	 
			Efficacy_CamelControls = unique(ModelRunsCopy$Efficacy_CamelControls)[EffCamelIndex]
			cat(paste0("Efficacy_CamelControls ", Efficacy_CamelControls, ":\t"))
			
			# subset model runs
			ModelRuns_subset_interim = ModelRunsCopy[	ModelRunsCopy$VacCampStrategy 			== "PROACTIVE" 				& 
							ModelRunsCopy$ImmunityDelay 			== 0 						& 
							ModelRunsCopy$Efficacy_CamelControls 	== Efficacy_CamelControls	, ]
			
			TimeIndex = 2
			for (TimeIndex in 1:length(unique(ModelRunsCopy$TimeSinceVaccination)))
			{
				TimeSinceVaccination = unique(ModelRunsCopy$TimeSinceVaccination)[TimeIndex]
				cat(paste0(TimeSinceVaccination, " years, "))
				
				# subset model runs
				ModelRuns_subset = ModelRuns_subset_interim[ModelRuns_subset_interim$TimeSinceVaccination == TimeSinceVaccination, ]
				
				if (any(c("FinalCaseDate_CF", "Peak") == ZVar_Chars[PlotNum])) LIMITS = NULL else LIMITS = c(0,1)
				
				PlotTitle = paste0("\n", TimesSinceVaccinationNames[TimeIndex], " wait")
				if (TimeIndex == 1)	PlotTitle = paste0("Camel control effectiveness = ", Efficacy_CamelControls, "%", PlotTitle)
				
				GPlot <- ggplot(ModelRuns_subset, aes_string(x = XVar_Char, y = YVar_Char, fill = paste0(ZVar_Chars[PlotNum], Stat))) + 
						geom_tile() + 
						ggtitle(PlotTitle) +
						labs(x = "", y = "") + 
						theme_minimal() + 
						theme(plot.title = element_text(size = 25), axis.text = element_text(size = 22.5), 
								legend.text = element_text(size = 25), legend.title = element_text(size = 29),
								legend.key.size = unit(2.25, 'cm'), legend.key.width = unit(2.2, 'cm')) + 
						guides(fill = FALSE) +  
						scale_fill_gradientn(
								colors = rev(matlab.like2(50)), # various colour candidates here: https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf #colors = heat.colors(50), #colors = hcl.colors(50, palette = "RdYlBu"), palette = "RdYlGn", palette = "YlOrRd"
								limits = LIMITS)
				
				assign(paste0("CamelEff_", EffCamelIndex, "_Time_", TimeIndex, "_Plot"),  GPlot	)	## use assign to use characters to individual geom_tile plots
				
				if (TimeIndex == length(unique(ModelRunsCopy$TimeSinceVaccination)))
				{
					GPlotWLegend <- GPlot + guides(fill = guide_legend(title = LegendNames[PlotNum]), limits = LIMITS) 
					legend <- get_legend(GPlotWLegend)
					assign(paste0("Legend_Row_", EffCamelIndex), legend)	## use assign to use characters to individual geom_tile plots
				}
			}
			cat(paste0("\n"))
		}
		cat(paste0("\n"))
		
		png(filename = file.path(here("Plots/"), paste0("Heatmaps_Proactive_", ZVar_Chars[PlotNum], Stat, ".png")), 
				res = 500, units = "in", width = 21, height = 15)
		
		LayoutMatrix = matrix(1:(length(unique(ModelRunsCopy$Efficacy_CamelControls)) * length(unique(ModelRunsCopy$TimeSinceVaccination))), 
				nrow = length(unique(ModelRunsCopy$Efficacy_CamelControls)), byrow = T)
		LayoutMatrix = cbind(LayoutMatrix, max(LayoutMatrix) + 1)
		
		grid.arrange(
				
				CamelEff_1_Time_1_Plot	,
				CamelEff_1_Time_2_Plot	,
				CamelEff_1_Time_3_Plot	,
				CamelEff_1_Time_4_Plot	,
				CamelEff_1_Time_5_Plot	,
				
				CamelEff_2_Time_1_Plot	,
				CamelEff_2_Time_2_Plot	,
				CamelEff_2_Time_3_Plot	,
				CamelEff_2_Time_4_Plot	,
				CamelEff_2_Time_5_Plot	,
				
				CamelEff_3_Time_1_Plot	,
				CamelEff_3_Time_2_Plot	,
				CamelEff_3_Time_3_Plot	,
				CamelEff_3_Time_4_Plot	,
				CamelEff_3_Time_5_Plot	,
				Legend_Row_3			,
				
				layout_matrix = LayoutMatrix,
				widths 	= c(rep(1, length(unique(ModelRunsCopy$TimeSinceVaccination))), 1.2),
				
				bottom 	= textGrob("Vaccine efficacy (%)"		, gp=gpar(fontsize = 35, font = 1), hjust = 0.85), 
				left 	= textGrob("Vaccine duration (years)"	, gp=gpar(fontsize = 35, font = 1), rot = 90)) ## would be better to use grid arrange with characters but leave for now.
		dev.off()
	}
}

CloseOpenPlotDevices()	



