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
head(ModelRuns, 200)
str(ModelRuns)
ListOfCounterfactuals = c("PropCasesAverted", "DeathsAverted_CF")

XVar_Char = "VaccineDuration"
YVar_Char = "TimeSinceVaccination"

ZVar_Chars 	= ListOfCounterfactuals
LegendNames	= c("Ratio of\ncases averted\n(proactive\nvs. reactive)", "Ratio of\ndeaths averted\n(proactive\nvs. reactive)")
str(ModelRuns$VaccineDuration)
str(ModelRuns)

# Make copy of model runs so can alter variables. E.g. non-uniform spacing between values in TimeSinceVaccination, or 0 being coded as no waning (i.e. infinite duration).
ModelRunsCopy = ModelRuns
ModelRunsCopy = ModelRunsCopy[ModelRunsCopy$Efficacy_Start 			%in% unique(ModelRunsCopy$Efficacy_Start)[c(2,10,18)], ]
ModelRunsCopy = ModelRunsCopy[ModelRunsCopy$Efficacy_CamelControls 	%in% unique(ModelRunsCopy$Efficacy_CamelControls)[c(1,2,4,6)], ]

### need ModelRunsCopy$VaccineDuration as factor
ModelRunsCopy$VaccineDuration[which(ModelRunsCopy$VaccineDuration == 0)] = Inf
ModelRunsCopy$Efficacy_Start 			= ModelRunsCopy$Efficacy_Start * 100
ModelRunsCopy$Efficacy_CamelControls 	= ModelRunsCopy$Efficacy_CamelControls * 100
ModelRunsCopy$VaccineDuration 			= as.factor(ModelRunsCopy$VaccineDuration		)
ModelRunsCopy$TimeSinceVaccination 		= as.factor(ModelRunsCopy$TimeSinceVaccination	)

unique(ModelRunsCopy$Efficacy_Start)
unique(ModelRunsCopy$VaccineDuration)

ImplementationDelay 	= 14
Efficacy_CamelControls 	= 0.3
ReactLevel 				= "HOSPITAL"
Efficacy_Start			= 50

unique(ModelRunsCopy$VaccineDuration)
unique(ModelRunsCopy$TimeSinceVaccination)
TimesSinceVaccinationNames = paste0(unique(ModelRunsCopy$TimeSinceVaccination), " year")

ContourValue 	= 1
ContourCol	 	= "black"

WindowsList = list(c(0, 548), c(0,180), c(181,364), c(365,548))
#WindowsList = list(c(0,180), c(181,364), c(365,548))
WaningInReactive = 1
for (WaningInReactive in 1)
{
	if (WaningInReactive == 0)	RW_String = "" else RW_String = "_RW"
	
	Window 		= WindowsList[[1]]
	WindowIndex = 1
	for (WindowIndex in 1:length(WindowsList))
	{
		cat(paste0("WindowIndex? ", WindowIndex, "\n"))
		Window	= WindowsList[[WindowIndex]]
		
		MinDay_Pruning = Window[1]
		MaxDay_Pruning = Window[2]
		
		if (MinDay_Pruning > 0 | MaxDay_Pruning < (MAX_ONSET_DAY - 1))	
			WindowString 	= paste0("_s", MinDay_Pruning, "_f", MaxDay_Pruning) else	WindowString 	= ""
		
		if (MinDay_Pruning == 0 	& MaxDay_Pruning == (MAX_ONSET_DAY - 1)	) LIMITS = c(0, 1.5) else 
		if (MinDay_Pruning == 0 	& MaxDay_Pruning == 180					) LIMITS = c(0, 13) else 
		if (MinDay_Pruning == 181 	& MaxDay_Pruning == 364					) LIMITS = c(0, 15) else 
		if (MinDay_Pruning == 365 	& MaxDay_Pruning == (MAX_ONSET_DAY - 1)	) LIMITS = c(0, 3)  
		

		for (ExpWaning in 1:0)
		{
			cat(paste0("ExpWaning? ", ExpWaning, "\n"))
			
			if (WindowIndex == 1) 
			{
				if (ExpWaning == 1) LIMITS = c(0, 1.5) else LIMITS = c(0, 2)
				
			} else 	if (WindowIndex == 2) LIMITS = c(0, 12.5) 	else
			if (WindowIndex == 3) LIMITS = c(0, 15) 	else
			if (WindowIndex == 4) LIMITS = c(0, 3) 		
			
			if (ExpWaning == 0)	WaningString 	= "_Hill" 						else WaningString 	= ""
			if (ExpWaning == 0)	DurString 		= "Vaccine half-life (years)" 	else DurString 		= "Vaccine duration of protection (years)"
			
			PlotNum = 1
			for (PlotNum in 1:length(ZVar_Chars))
			{
				cat(paste0("PlotNum ", PlotNum, ": ", ZVar_Chars[PlotNum], "\n"))
				
				EffCamelIndex = 1
				for (EffCamelIndex in 1:length(unique(ModelRunsCopy$Efficacy_CamelControls)))
				{	 
					Efficacy_CamelControls = unique(ModelRunsCopy$Efficacy_CamelControls)[EffCamelIndex]
					cat(paste0("Efficacy_CamelControls ", Efficacy_CamelControls, ":\t"))
					
					ImplementationDelay = 28
					for (ImplementationDelay in c(8, 28))
					{
						cat(paste0("ImplementationDelay ", ImplementationDelay, ":\t"))
						
						Stat = "_Mean"
						for (Stat in c("_Mean"))
						{
#							for (ReactLevelIndex in 1:length(unique(ModelRunsCopy$ReactLevel)))
							{
								ReactLevel = unique(ModelRunsCopy$ReactLevel)[ReactLevelIndex]
								cat(paste0("ReactLevel ", ReactLevel, ":\t"))
								
								VariableInQuestion 	= paste0(ZVar_Chars[PlotNum], Stat)
								
								EffHumanIndex = 3
								for (EffHumanIndex in 1:length(unique(ModelRunsCopy$Efficacy_Start)))
								{
									Efficacy_Start = unique(ModelRunsCopy$Efficacy_Start)[EffHumanIndex]
									
									# Make reactive subset
									ModelRuns_Reactive = ModelRunsCopy[ModelRunsCopy$VacCampStrategy 	== "REACTIVE" 		& 
													ModelRunsCopy$ReactLevel 				== ReactLevel					&
													ModelRunsCopy$WaningInReactive			== WaningInReactive				&
													ModelRunsCopy$ImmunityDelay 			== 14 							& 
													ModelRunsCopy$ImplementationDelay 		== ImplementationDelay			& 
													ModelRunsCopy$MinDay_Pruning 			== MinDay_Pruning 				& 
													ModelRunsCopy$MaxDay_Pruning 			== MaxDay_Pruning 				& 
													ModelRunsCopy$Efficacy_Start 			== Efficacy_Start 				& 
													ModelRunsCopy$Efficacy_CamelControls 	== Efficacy_CamelControls, 		]
									
									head(ModelRuns_Reactive, 100)
									
									if (WaningInReactive == 1)
										ModelRuns_Reactive = ModelRuns_Reactive[ModelRuns_Reactive$VaccineDuration != Inf 	&
														ModelRuns_Reactive$ExpWaning 				== ExpWaning ,]
									
									# Now make proactive subset and subtract reactive subset.
									ModelRuns_Proactive = ModelRunsCopy[ModelRunsCopy$VacCampStrategy 	== "PROACTIVE" 		& 
													ModelRunsCopy$ImmunityDelay 			== 0 							& 
													ModelRunsCopy$VaccineDuration 			!= Inf 							& 
													ModelRunsCopy$Efficacy_Start 			== Efficacy_Start 				& 
													ModelRunsCopy$ExpWaning 				== ExpWaning 				& 
													ModelRunsCopy$MinDay_Pruning 			== MinDay_Pruning 			& 
													ModelRunsCopy$MaxDay_Pruning 			== MaxDay_Pruning 			& 
													ModelRunsCopy$Efficacy_CamelControls 	== Efficacy_CamelControls, 		]
									
									cbind(	ModelRuns_Proactive[,c("Efficacy_Start", "VaccineDuration", "TimeSinceVaccination")] , 
											ModelRuns_Reactive[,c("Efficacy_Start", "VaccineDuration", "ImplementationDelay","ReactLevel")], Pro = ModelRuns_Proactive[,VariableInQuestion] , Re = ModelRuns_Reactive[,VariableInQuestion])
									length(  ModelRuns_Reactive[,VariableInQuestion])
									length(  ModelRuns_Proactive[,VariableInQuestion])
									
									
									# divide proactive by reactive (for VariableInQuestion) to obtain ratio.
									ModelRuns_Proactive[,VariableInQuestion] = ModelRuns_Proactive[,VariableInQuestion] / ModelRuns_Reactive[,VariableInQuestion]
									max(ModelRuns_Proactive[,VariableInQuestion])	## proactive / reactive
									max(1/ModelRuns_Proactive[,VariableInQuestion])	## reactive / proactive
									
									PlotTitle = paste0("\n", "Vaccine efficacy = ", Efficacy_Start, "%")
									if (EffHumanIndex == 1)	PlotTitle = paste0(ReactLevel, " LEVEL", PlotTitle)
									
									GPlot <- ggplot(ModelRuns_Proactive, 
													aes_string(x = XVar_Char, y = YVar_Char, z = VariableInQuestion, fill = VariableInQuestion)) + 
											geom_tile() + 
											ggtitle(PlotTitle) +
											labs(x = "", y = "") + 
											theme_minimal() +  
											theme(plot.title = element_text(size = 30), axis.text = element_text(size = 25), 
													legend.text = element_text(size = 25), legend.title = element_text(size = 31),
													legend.key.size = unit(2.75, 'cm'), legend.key.width = unit(2.2, 'cm'), 
													panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
											guides(fill = FALSE) +  
											scale_fill_gradientn(
													#colors = rev(matlab.like2(50)), # various colour candidates here: https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf #colors = heat.colors(50), #colors = hcl.colors(50, palette = "RdYlBu"), palette = "RdYlGn", palette = "YlOrRd"
													colors = matlab.like2(50), # various colour candidates here: https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf #colors = heat.colors(50), #colors = hcl.colors(50, palette = "RdYlBu"), palette = "RdYlGn", palette = "YlOrRd"
													limits = LIMITS, breaks = seq(LIMITS[1], LIMITS[2], length.out = 11))  #+ stat_contour()  
									
									XValueIndex = 2
									for (XValueIndex in 1:length(unique(ModelRuns_Proactive[,XVar_Char])))
									{
										XValue_Fact = unique(ModelRuns_Proactive[,XVar_Char])[XValueIndex]
										XValue_Num  = unique(as.numeric(ModelRuns_Proactive[,XVar_Char]))[XValueIndex]
										if (any(ModelRuns_Proactive[ModelRuns_Proactive[,XVar_Char] == XValue_Fact, VariableInQuestion] > ContourValue))
											if (any(ModelRuns_Proactive[ModelRuns_Proactive[,XVar_Char] == XValue_Fact, VariableInQuestion] < ContourValue)) ## need this second condition otherwise contour doesn't quite work
											{
												YValue_Num 	= which.min(ModelRuns_Proactive[ModelRuns_Proactive[,XVar_Char] == XValue_Fact, VariableInQuestion] > ContourValue)
												GPlot <- GPlot + geom_segment(x = XValue_Num - 0.5, y = YValue_Num - 0.5, xend = XValue_Num + 0.5, yend = YValue_Num - 0.5, 
														size = 3, colour = ContourCol)
											}
									}
									YValueIndex = 1
									for (YValueIndex in 1:length(unique(ModelRuns_Proactive[,YVar_Char])))
									{
										YValue_Fact = unique(ModelRuns_Proactive[,YVar_Char])[YValueIndex]
										YValue_Num  = unique(as.numeric(ModelRuns_Proactive[,YVar_Char]))[YValueIndex]
										if (any(ModelRuns_Proactive[ModelRuns_Proactive[,YVar_Char] == YValue_Fact, VariableInQuestion] > ContourValue))
											if (any(ModelRuns_Proactive[ModelRuns_Proactive[,YVar_Char] == YValue_Fact, VariableInQuestion] < ContourValue)) ## need this second condition otherwise contour doesn't quite work
											{
												XValueIndex 	= which.min(ModelRuns_Proactive[ModelRuns_Proactive[,YVar_Char] == YValue_Fact, VariableInQuestion] > ContourValue)
												XValue_Fact 	= unique(ModelRuns_Proactive[,XVar_Char])[XValueIndex]
												XValue_Num  	= unique(as.numeric(ModelRuns_Proactive[,XVar_Char]))[XValueIndex]
												GPlot <- GPlot + geom_segment(x = XValue_Num + 0.5, y = YValue_Num - 0.5, xend = XValue_Num + 0.5, yend = YValue_Num + 0.5, 
														size = 3, colour = ContourCol)
											}
									}
									
									assign(paste0("ReactLevel_", ReactLevelIndex, "_VE_", EffHumanIndex, "_Plot"),  GPlot	)	## use assign to use characters to individual geom_tile plots
									
									if (EffHumanIndex == length(unique(ModelRunsCopy$Efficacy_Start)))
									{
										GPlotWLegend <- GPlot + guides(fill=guide_legend(title = LegendNames[PlotNum]), limits = LIMITS) 
										legend <- get_legend(GPlotWLegend)
										assign(paste0("Legend_Row_", ReactLevelIndex), legend)	## use assign to use characters to individual geom_tile plots
									}
								}
								cat(paste0("\n"))
							}
							cat(paste0("\n"))
							
							LayoutMatrix = matrix(1:(length(unique(ModelRunsCopy$Efficacy_Start)) * length(unique(ModelRunsCopy$ReactLevel))), 
									nrow = length(unique(ModelRunsCopy$ReactLevel)), byrow = T)
							LayoutMatrix = cbind(LayoutMatrix, max(LayoutMatrix) + 1)
							
							Filename = paste0("ProactiveVsReactive_CamelEff_", Efficacy_CamelControls, "_Delay_", ImplementationDelay, "_", ZVar_Chars[PlotNum], Stat, WindowString, WaningString, RW_String)
							png(filename = file.path(here("Plots/"), paste0(Filename, ".png")),	res = 500, units = "in", width = 21, height = 17)
#							pdf(file = file.path(here(), paste0(Filename, ".pdf")),	width = 21, height = 17)
							
							grid.arrange(
									
									ReactLevel_1_VE_1_Plot	,
									ReactLevel_1_VE_2_Plot	,
									ReactLevel_1_VE_3_Plot	,
									
									ReactLevel_2_VE_1_Plot	,
									ReactLevel_2_VE_2_Plot	,
									ReactLevel_2_VE_3_Plot	,
									
									ReactLevel_3_VE_1_Plot	,
									ReactLevel_3_VE_2_Plot	,
									ReactLevel_3_VE_3_Plot	,
									Legend_Row_1			,
									
									layout_matrix = LayoutMatrix,
									widths 	= c(rep(1, length(unique(ModelRunsCopy$ReactLevel))), 0.8),
									
									nrow 	= length(unique(ModelRunsCopy$ReactLevel)), 
									#top 	= textGrob(paste0("Ratio of cases averted (Proactive vs. reactive)"), gp = gpar(fontsize = 45, font = 1)), 
									bottom 	= textGrob(DurString								, gp=gpar(fontsize = 35, font = 1), hjust = 0.7), 
									left 	= textGrob("Vaccination-to-outbreak wait (years)"	, gp=gpar(fontsize = 35, font = 1), rot = 90)) ## would be better to use grid arrange with characters but leave for now.
							dev.off()
						}
					}
				}
			}
		}
	}
}

CloseOpenPlotDevices()

