
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

XVar_Char = "ImplementationDelay"
YVar_Char = "Efficacy_Start"

ReactLevel 	= "HOSPITAL"; 
ZVar_Chars 	= c("PropCasesAverted"	, "DeathsAverted_CF")
PlotNames 	= c("cases averted"		, "deaths averted"	)
LegendNames	= c("Proportion of\ncases averted", "Proportion of\ndeaths averted")
#ZVar_Chars 	= c(ListOfCounterfactuals, "Peak_Reduction")
#LegendNames	= c("Proportion of\ncases averted"	, "Proportion of\ncases averted\n(healthcare\nworkers only)" , "Proportion of\ncases averted\n(excluding\nhealthcare workers)", 
#		"Proportion of\ndeaths averted"	, "Proportion of\ndeaths averted\n(healthcare\nworkers only)", "Proportion of\ndeaths averted\n(excluding\nhealthcare workers)", 
#		"Final\ncase day"				, "Epidemic\npeak size", "Proportional\nreduction in\nepidemic peak")
Stat 		= "_Mean"

CamelEffsSubset = unique(ModelRuns$Efficacy_CamelControls)[c(1,4,6)]

PlotNum = 1
for (PlotNum in 1:length(ZVar_Chars))
{
	if (ZVar_Chars[PlotNum] == "FinalCaseDate_CF") LIMITS = NULL else LIMITS = c(0,1)
	
	for (Stat in c("_Mean", "_LowerCrI", "_UpperCrI"))
	{
		for (EffCamelIndex in 1:length(CamelEffsSubset))
		{	 
			Efficacy_CamelControls = CamelEffsSubset[EffCamelIndex]
			for (ReactLevel in unique(ModelRuns$ReactLevel))
			{
				# subset model runs
				ModelRuns_subset = ModelRuns[ModelRuns$VacCampStrategy == "REACTIVE" & 
								ModelRuns$ReactLevel == ReactLevel & ModelRuns$Efficacy_CamelControls == Efficacy_CamelControls & 
								ModelRuns$Efficacy_Start != 0 & 
								ModelRuns$ImmunityDelay == 14, ]
				ModelRuns_subset$Efficacy_Start = ModelRuns_subset$Efficacy_Start * 100
				
				cat(paste0("ReactLevel ", ReactLevel, ", "))
				
				if (ReactLevel == "HOSPITAL")	ReactLevel_String = "Hospital"
				if (ReactLevel == "REGIONAL")	ReactLevel_String = "Regional"
				if (ReactLevel == "NATIONAL")	ReactLevel_String = "National"
				
				PlotTitle = paste0("\n", ReactLevel_String, " level")
				if (ReactLevel == "HOSPITAL")	PlotTitle = paste0("Camel control effectiveness = ", Efficacy_CamelControls * 100, "%", PlotTitle)
				
				# make heatmap
				GPlot <- ggplot(ModelRuns_subset, aes_string(x = XVar_Char, y = YVar_Char, fill = paste0(ZVar_Chars[PlotNum], Stat))) + 
						geom_tile() + 
						ggtitle(PlotTitle) + 
						labs(y = "", x = "") + 
						theme_minimal() + 
						theme(plot.title = element_text(size = 30), axis.text = element_text(size = 25), 
								legend.text = element_text(size = 25), legend.title = element_text(size = 31),
								legend.key.size = unit(2.25, 'cm'), legend.key.width = unit(2.2, 'cm')) + 
						guides(fill = FALSE) + 
						scale_fill_gradientn(
								colors = rev(matlab.like2(50)), # various colour candidates here: https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf #colors = heat.colors(50), #colors = hcl.colors(50, palette = "RdYlBu"), palette = "RdYlGn", palette = "YlOrRd"
								limits = LIMITS)
				
				# assign it to nice name
				assign(paste0("CamelEff_", EffCamelIndex, "_", ReactLevel, "_Plot"), GPlot)
				
				if (ReactLevel == "NATIONAL")
				{
					GPlotWLegend <- GPlot + guides(fill = guide_legend(title = LegendNames[PlotNum]), limits = LIMITS) 
					legend <- get_legend(GPlotWLegend)
					assign(paste0("Legend_Row_", EffCamelIndex),  legend)	
				}
			}
			cat(paste0(Stat, ", "))
		}
		
		png(filename = file.path(here("Plots/"), paste0("Heatmap_Reactive", Stat, "_", ZVar_Chars[PlotNum], ".png")), 
				res = 500, units = "in", width = 21, height = 17)
		
		LayoutMatrix = matrix(1:(length(CamelEffsSubset) * length(unique(ModelRuns$ReactLevel))), 
				nrow = length(CamelEffsSubset), byrow = T)
		LayoutMatrix = cbind(LayoutMatrix, max(LayoutMatrix) + 1)
		
		grid.arrange(
				
				CamelEff_1_HOSPITAL_Plot	, 
				CamelEff_1_REGIONAL_Plot	, 
				CamelEff_1_NATIONAL_Plot	, 
				CamelEff_2_HOSPITAL_Plot	, 
				CamelEff_2_REGIONAL_Plot	, 
				CamelEff_2_NATIONAL_Plot	, 
				CamelEff_3_HOSPITAL_Plot	, 
				CamelEff_3_REGIONAL_Plot	, 
				CamelEff_3_NATIONAL_Plot	, 
				Legend_Row_1					,
				
				layout_matrix = LayoutMatrix,
				widths 	= c(1,1,1,0.8)			,
				
				nrow 	= length(CamelEffsSubset), 
				bottom 	= textGrob("React time (days)"		, gp = gpar(fontsize = 35, font = 1), hjust = 0.97), 
				left 	= textGrob("Vaccine efficacy (%)"	, gp = gpar(fontsize = 35, font = 1), rot = 90)) ## would be better to use grid arrange with characters but leave for now.
		dev.off()
		
		cat(paste0("\n"))
	}
}
class(GPlot)

CloseOpenPlotDevices()	
