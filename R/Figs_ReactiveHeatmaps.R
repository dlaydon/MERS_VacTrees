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
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
PlotsDirectory		= file.path(ProjectDirectory, "Plots"		)

source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))
source(file.path(R_ScriptDirectory, "LayoutPlottingFunctions.R"))

# Import summary model runs
ModelRuns = read.table(file = file.path(ProjectDirectory, paste0("ModelRunsSummary.txt")), header = T, sep = "\t")
ModelRuns = ModelRuns[ModelRuns$DoTriggers == 0, ]
str(ModelRuns)

ZVar_Chars 	= c("PropCasesAverted"	, "DeathsAverted_CF", "Cases_DirectlyProtected"	, "Cases_IndirectlyProtected", "Deaths_DirectlyProtected", "Deaths_IndirectlyProtected")
LegendNames	= c("Proportion of\ncases averted", "Proportion of\ndeaths averted", "Proportion of\ncases averted\n(directly)", "Proportion of\ncases averted\n(indirectly)", "Proportion of\ndeaths averted\n(directly)", "Proportion of\ndeaths averted\n(indirectly)")

CamelEffsSubset = unique(ModelRuns$Efficacy_CamelControls)[c(1,4,6)]
DurationSubset 	= unique(ModelRuns$VaccineDuration)[c(7,6,5,4,2)]

WindowsList = list(c(0, 548), c(0,180), c(181,364), c(365,548))

RW_String 	= "_RW"
WaningInReactive = 1

XVar_Char = "ImplementationDelay"
YVar_Char = "Efficacy_Start"

MinDay_Pruning = Window[1]
MaxDay_Pruning = Window[2]
if (MinDay_Pruning > 0 | MaxDay_Pruning < (MAX_ONSET_DAY - 1))	
	WindowString = paste0("_s", MinDay_Pruning, "_f", MaxDay_Pruning) else	WindowString = "" 
WindowString = ""

EffCamelIndex = 1
for (EffCamelIndex in 1:length(CamelEffsSubset))
{	 
	Efficacy_CamelControls = CamelEffsSubset[EffCamelIndex]
	cat(paste0("Efficacy_CamelControls ", Efficacy_CamelControls, ", \t"))
	
	ExpWaning = 0
	for (ExpWaning in 1:0)
	{
		if (ExpWaning == 0)	HillString 	= "_Hill" else HillString = ""
		
		cat(paste0("ExpWaning ", ExpWaning, ", "))
		
		PlotNum = 2
		for (PlotNum in 1:length(ZVar_Chars))
		{
			if (ZVar_Chars[PlotNum] == "FinalCaseDate_CF") LIMITS = NULL else LIMITS = c(0,1)
			
			cat(paste0(ZVar_Chars[PlotNum], " "))
			
			Stat 		= "_Mean"
			for (Stat in c("_Mean", "_LowerCrI", "_UpperCrI"))
			{
				cat(paste0(Stat, "\n"))
				Duration = DurationSubset[1]
				for (DurationIndex in 1:length(DurationSubset))
				{
					Duration = DurationSubset[DurationIndex]
					cat(paste0("Duration ", Duration, ", \t "))
					
					ReactLevel 	= "HOSPITAL"; 
					ReactLevel 	= "REGIONAL"; 
					ReactLevel 	= "NATIONAL"; 
					for (ReactLevel in unique(ModelRuns$ReactLevel))
					{
						cat(paste0(ReactLevel, ", \t"))
						
						ModelRuns_subset = ModelRuns[ModelRuns$VacCampStrategy 	== "REACTIVE" 				& 
										ModelRuns$ReactLevel 					== ReactLevel 				& 
										ModelRuns$WaningInReactive				== WaningInReactive			&
										ModelRuns$Efficacy_CamelControls 		== Efficacy_CamelControls 	& 
										ModelRuns$Efficacy_Start 				!= 0 						& 
										ModelRuns$MinDay_Pruning 				== MinDay_Pruning 			& 
										ModelRuns$MaxDay_Pruning 				== MaxDay_Pruning 			& 
										ModelRuns$ImmunityDelay 				== 14						, ]
						
						if (WaningInReactive == 1)
							ModelRuns_subset = ModelRuns_subset[ModelRuns_subset$VaccineDuration == Duration 	&
											ModelRuns_subset$ExpWaning 				== ExpWaning ,]					
						
						
						ReactLevel_String = str_to_sentence(ReactLevel, locale = "en")
						PlotTitle = paste0("\n", ReactLevel_String, " level")
						
						if (ExpWaning == 0			)	DurString 	= "Half-life" 	else DurString = "Mean duration of protection"
						if (Duration == 1			)	YearOrYears = "year" 		else YearOrYears = "years"
						if (ReactLevel == "HOSPITAL")	PlotTitle 	= paste0(DurString, " = ", Duration, " ", YearOrYears, PlotTitle)
						
						# make heatmap
						GPlot <- ggplot(ModelRuns_subset, aes_string(x = XVar_Char, y = YVar_Char, fill = paste0(ZVar_Chars[PlotNum], Stat))) + 
								geom_tile() + 
								ggtitle(PlotTitle) + 
								labs(y = "", x = "") + 
								theme_minimal() + 
								theme(plot.title = element_text(size = 30), axis.text = element_text(size = 25), 
										legend.text = element_text(size = 25), legend.title = element_text(size = 31),
										legend.key.size = unit(2.25, 'cm'), legend.key.width = unit(2.2, 'cm'), 
										panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
								guides(fill = "none") + 
								scale_fill_gradientn(colors = rev(matlab.like2(50)), limits = LIMITS, breaks = seq(LIMITS[1], LIMITS[2], by = 0.1)) # various colour candidates here: https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf #colors = heat.colors(50), #colors = hcl.colors(50, palette = "RdYlBu"), palette = "RdYlGn", palette = "YlOrRd"
						
						# assign it to nice name
						assign(paste0("Dur_", Duration, "_", ReactLevel, "_Plot"), GPlot)
						
						if (ReactLevel == "NATIONAL")
						{
							GPlotWLegend <- GPlot + guides(fill = guide_legend(title = LegendNames[PlotNum]), limits = LIMITS, breaks = seq(LIMITS[1], LIMITS[2], by = 0.1)) 
							legend <- get_legend(GPlotWLegend)
							assign(paste0("Legend_Row_", EffCamelIndex),  legend)	
						}
					}
					cat(paste0("\n"))
				}
				
				Filename = paste0("ReactiveW_", ZVar_Chars[PlotNum], WindowString, HillString, "CamelEff_", Efficacy_CamelControls, Stat)
				png(filename 	= file.path(here("Plots/"), paste0(Filename, ".png")),	res = 500, units = "in"	, width = 21, height = 27.5)
				#pdf(file 		= file.path(here("Plots/"), paste0(Filename, ".pdf"))							, width = 21, height = 17)
				
				LayoutMatrix = matrix(1:(length(DurationSubset) * length(unique(ModelRuns$ReactLevel))), nrow = length(DurationSubset), byrow = T)
				LayoutMatrix = cbind(LayoutMatrix, max(LayoutMatrix) + 1)
				
				grid.arrange(
						
						Dur_1_HOSPITAL_Plot	, 
						Dur_1_REGIONAL_Plot	, 
						Dur_1_NATIONAL_Plot	, 
						
						Dur_2_HOSPITAL_Plot	, 
						Dur_2_REGIONAL_Plot	, 
						Dur_2_NATIONAL_Plot	, 
						
						Dur_5_HOSPITAL_Plot	, 
						Dur_5_REGIONAL_Plot	, 
						Dur_5_NATIONAL_Plot	, 
						
						Dur_10_HOSPITAL_Plot	, 
						Dur_10_REGIONAL_Plot	, 
						Dur_10_NATIONAL_Plot	, 
						
						Dur_20_HOSPITAL_Plot	, 
						Dur_20_REGIONAL_Plot	, 
						Dur_20_NATIONAL_Plot	, 
						Legend_Row_1			,
						
						layout_matrix = LayoutMatrix,
						widths 	= c(1,1,1,0.8)			,
						
						nrow 	= length(DurationSubset), 
						bottom 	= textGrob("Reaction time (days)"	, gp = gpar(fontsize = 35, font = 1), hjust = 0.875), 
						left 	= textGrob("Vaccine efficacy (%)"	, gp = gpar(fontsize = 35, font = 1), rot = 90)) ## would be better to use grid arrange with characters but leave for now.
				dev.off()
			}
		}
	}	
}

CloseOpenPlotDevices()
warnings()

