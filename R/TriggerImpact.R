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
OrigMAI = par ("mai") ### good to record if using layout functions
OrigMAR = par ("mar") ### good to record if using layout functions

ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppRootDirectory 	= file.path(ProjectDirectory, "MERS_Vac"	, "MERS_Vac") 
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
RawDataDirectory 	= file.path(ProjectDirectory, "Data"		)
source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))

# Import processed data for Cpp
DATA 	= read.table(file = file.path(RawDataDirectory, "MERS_forCpp.txt"), header = TRUE)
MR = read.table(file = file.path(ProjectDirectory, "ModelRunsSummary.txt"), header = T, sep = "\t")
MR = MR[MR$VacCampStrategy == "REACTIVE",]
MR = MR[MR$ExpWaning == 1,]
MR = MR[MR$Efficacy_CamelControls == 0,]
MR = MR[MR$WaningInReactive == 1,]
MR = MR[MR$MinDay_Pruning == 0,]
MR = MR[MR$MaxDay_Pruning == 548,]
MR = MR[MR$ImmunityDelay == 14,]
MR = MR[MR$ImplementationDelay %in% c(14, 28),]

# rename (make copy) with nicer variable names
MR$Level 		= MR$ReactLevel
MR$Efficacy 	= MR$Efficacy_Start
MR$Duration 	= MR$VaccineDuration
MR$Threshold 	= MR$TrigThreshold
MR$Timeframe 	= MR$TrigTimeframe
MR$ReactionTime = MR$ImplementationDelay
MR$Impact 		= MR$PropCasesAverted_Mean

# select relevant runs and variables
MRNoTrig 	= MR[MR$DoTriggers == 0, c("Level", "Efficacy", "Duration", "Threshold", "Timeframe", "ReactionTime", "Impact")]
MR 			= MR[MR$DoTriggers == 1, c("Level", "Efficacy", "Duration", "Threshold", "Timeframe", "ReactionTime", "Impact")]
MR$Ratio 	= NA

# populate Ratio
for (Level in unique(MR$Level))
	for (Efficacy in unique(MR$Efficacy))
		for (Duration in unique(MR$Duration))
			for (Threshold in unique(MR$Threshold))
				for (Timeframe in unique(MR$Timeframe))
					for (ReactionTime in unique(MR$ReactionTime))
{
	IndexInTrigRuns = which(MR$Level == Level & MR$Efficacy == Efficacy & MR$Duration == Duration & 
					MR$Threshold == Threshold & MR$Timeframe == Timeframe & MR$ReactionTime == ReactionTime)
	
	IndexInNonTrigRuns = which(MRNoTrig$Level == Level & MRNoTrig$Efficacy == Efficacy & MRNoTrig$Duration == Duration & 
					MRNoTrig$ReactionTime == ReactionTime)
	
	MR$Ratio[IndexInTrigRuns] = MR$Impact[IndexInTrigRuns] / MRNoTrig$Impact[IndexInNonTrigRuns]
}
						
str(MRNoTrig)
str(MR)


Efficacy = 0.9
Level = "HOSPITAL"
Threshold = 5
Timeframe = 14

Duration = 10
ReactionTime = 28
LIMITS = c(0,1)
EffSub = c(0.5, 0.9)

for (ReactionTime in unique(MR$ReactionTime))
	for (EfficacyIndex in 1:length(EffSub))
		for (Level in unique(MR$Level))
		{
			Efficacy = EffSub[EfficacyIndex]
			
			MR_Sub = MR[MR$Level == Level & MR$Efficacy == Efficacy & MR$Duration == Duration & MR$ReactionTime == ReactionTime, ]
			MR_Sub$Threshold = as.factor(MR_Sub$Threshold)
			MR_Sub$Timeframe = as.factor(MR_Sub$Timeframe)
			
			ReactLevel_String = str_to_sentence(Level, locale = "en")
			PlotTitle = paste0("\n", ReactLevel_String, " level")
			
			if (Level == "HOSPITAL") PlotTitle 	= paste0("Efficacy = ", Efficacy * 100, "%, ", ReactionTime, "-day reaction time", PlotTitle)
			
			GPlot <- ggplot(MR_Sub, aes(x = Threshold, y = Timeframe, fill = Ratio)) + 
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
			assign(paste0("Eff_", EfficacyIndex, "_", Level, "_ReactTime_", ReactionTime, "_Plot"), GPlot)
			
			if (Level == "NATIONAL")
			{
				GPlotWLegend <- GPlot + guides(fill = guide_legend(title = "Ratio"), limits = LIMITS, breaks = seq(LIMITS[1], LIMITS[2], by = 0.1)) 
				legend <- get_legend(GPlotWLegend)
				assign(paste0("Legend_Row_", 1),  legend)	
			}
			GPlot
		}


png(filename = here("Plots/TrigRatios.png"), res = 500, units = "in", width = 21, height = 23)

LayoutMatrix = matrix(1:(length(EffSub) * length(unique(MR$Level)) * length(unique(MR$ReactionTime))), ncol = length(unique(MR$Level)), byrow = T)
LayoutMatrix = cbind(LayoutMatrix, max(LayoutMatrix) + 1)

grid.arrange(
		
		Eff_1_HOSPITAL_ReactTime_14_Plot	, 
		Eff_1_REGIONAL_ReactTime_14_Plot	, 
		Eff_1_NATIONAL_ReactTime_14_Plot	, 
		
		Eff_1_HOSPITAL_ReactTime_28_Plot	, 
		Eff_1_REGIONAL_ReactTime_28_Plot	, 
		Eff_1_NATIONAL_ReactTime_28_Plot	, 
		
		Eff_2_HOSPITAL_ReactTime_14_Plot	, 
		Eff_2_REGIONAL_ReactTime_14_Plot	, 
		Eff_2_NATIONAL_ReactTime_14_Plot	, 
		
		Eff_2_HOSPITAL_ReactTime_28_Plot	, 
		Eff_2_REGIONAL_ReactTime_28_Plot	, 
		Eff_2_NATIONAL_ReactTime_28_Plot	, 
		
		Legend_Row_1			,
		
		layout_matrix = LayoutMatrix,
		widths 	= c(1,1,1,0.7)			,
		
		bottom 	= textGrob("Threshold (number of infections)"	, gp = gpar(fontsize = 35, font = 1), hjust = 0.70), 
		left 	= textGrob("Timeframe (days)"					, gp = gpar(fontsize = 35, font = 1), rot = 90)) ## would be better to use grid arrange with characters but leave for now.
dev.off()


