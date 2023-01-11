## LOAD PACKAGES ##
#install.packages(c("igraph", "adegenet", "RColorBrewer"))
library(igraph)
library(adegenet)
library(RColorBrewer)
require(here)
library(grid)
library(gridExtra)

## Want 2x5 figure showing epidemic curves and multiple who-infected-whom trees.
# 1st row - epidemic curves of cases, then deaths, with cases and deaths among HCWs highlighted. 
# 2nd row:
	# 1st column: no vaccine; 
	# 2nd column: regional reactive with 60% efficacy, 20-year mean duration (exponential waning), 14-day reaction time
	# 3rd column: proactive 60% efficacy, 20-year mean duration (exponential waning), 6-month wait; 
	# 4th column: Camels vaccinated (say Eff = 50%)

## Add #defines from Cpp code
AnimalReservoir		= -1 # note defined as 0 in Cpp but -1 here
SameHosptial		= 0 # i.e. transmission within same hospital
SameRegion			= 1 # i.e. transmission within same region but between different hospitals
DifferentRegion		= 2 # i.e. transmission between regions

# LOAD DATA ##
ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppRootDirectory 	= file.path(ProjectDirectory, "MERS_Vac"	, "MERS_Vac") 
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		)
RawDataDirectory 	= file.path(ProjectDirectory, "Data"		)
PlotsDirectory		= file.path(ProjectDirectory, "Plots"		)

source(here("R/DirectoryFunctions.R"))
source(here("R/LayoutPlottingFunctions.R"))

DATA 		= read.table(file = file.path(RawDataDirectory, "MERS_forCpp.txt"), header = TRUE)
Day_0  		= as.Date		("2013-01-01") #### 1st Jan 2013 is Day zero. 
Day_Final  	= Day_0 + max(DATA$onset) # note Different from Day_Final in MakeLineListForCpp.R script.
LastDay_Data_Int 	= max(DATA$onset)

## create vectors of for each strata you care about

Onset_Week = floor(DATA$onset / 7)

range(DATA$onset)
range(DATA$onset) + Day_0

AllCases_Inc 	= rep(NA, max(Onset_Week)); 
AllDeaths_Inc 	= rep(NA, max(Onset_Week)); 
Cases_HCW 		= rep(NA, max(Onset_Week)); 
Cases_nHCW 		= rep(NA, max(Onset_Week)); 
Deaths_HCW 		= rep(NA, max(Onset_Week)); 
Deaths_nHCW 	= rep(NA, max(Onset_Week)); 

table(Onset_Week)
table(Onset_Week[which(DATA$Dead == 1)])

AllCases_Inc	[as.numeric(names(table(Onset_Week)))] 							= table(Onset_Week)
AllDeaths_Inc	[as.numeric(names(table(Onset_Week[which(DATA$Dead == 1)])))] 	= table(Onset_Week[which(DATA$Dead == 1)])
Cases_HCW		[as.numeric(names(table(Onset_Week[which(DATA$HCW == 1)])))] 	= table(Onset_Week[which(DATA$HCW == 1)])
Cases_nHCW		[as.numeric(names(table(Onset_Week[which(DATA$HCW == 0)])))] 	= table(Onset_Week[which(DATA$HCW == 0)])
Deaths_HCW		[as.numeric(names(table(Onset_Week[which(DATA$Dead == 1 & DATA$HCW == 1)])))] 	= table(Onset_Week[which(DATA$Dead == 1 & DATA$HCW == 1)])
Deaths_nHCW		[as.numeric(names(table(Onset_Week[which(DATA$Dead == 1 & DATA$HCW == 0)])))] 	= table(Onset_Week[which(DATA$Dead == 1 & DATA$HCW == 0)])

length(which(DATA$Dead == 1 & DATA$HCW == 1))
length(which(DATA$Dead == 1 & DATA$HCW == 0))

length(which(DATA$HCW == 0))

DateLabels 			= c("Jan '13", "Apr '13", "Jul '13", "Oct '13", "Jan '14", "Apr '14", "Jul '14")
DatesForAxis 		= as.Date(c("2013-01-01", "2013-04-01", "2013-07-01", "2013-10-01", "2014-01-01", "2014-04-01", "2014-07-01"))
DatesForAxis_Int 	= round(as.numeric(DatesForAxis - Day_0) / 7)
DatesForAxis_Int[1] = 1 # not zero 

# import counterfactual trees for 2nd row of figure

# 1st column: no vaccine; 
Trees_noCFs 	<- read.table(here("Plots/ExampleTrees_nonCF.txt"), header = TRUE)
str(Trees_noCFs)
head(Trees_noCFs)

# 2nd column: regional reactive with 60% efficacy, 20-year mean duration (exponential waning), 14-day reaction time
ModelRun_Re_ModEff_SlowReact = DefineModelRuns(IncludeCompletedRunsOnly = TRUE, WaningInReactiveAndNot = 1, VaccineDurations = 20, 	
		VacCampStrategies = "REACTIVE", ReactLevels = "REGIONAL", Efficacies_Start = 0.6, ImplementationDelays = 14, ImmunityDelays = 14) 
Trees_Re_ModEff_SlowReact 	= read.table(file = file.path(CppOutputDirectory, paste0("Trees_CF_", ModelRun_Re_ModEff_SlowReact$OutputFolderNames, ".txt")), header = TRUE)

# 3rd column: proactive 60% efficacy, 20-year mean duration (exponential waning), 6-month wait; 
ModelRun_Pro_ModEff_ShortLag 	= DefineModelRuns(IncludeCompletedRunsOnly = TRUE, VacCampStrategies = "PROACTIVE", Efficacies_Start = 0.6, VaccineDurations = 5, TimesSinceVaccination = 0.5) 
Trees_Pro_ModEff_ShortLag 		= read.table(file = file.path(CppOutputDirectory, paste0("Trees_CF_", ModelRun_Pro_ModEff_ShortLag$OutputFolderNames, ".txt")), header = TRUE)

# 4th column: Camels vaccinated (say Eff = 50%)
ModelRun_CamelsOnly = DefineModelRuns(IncludeCompletedRunsOnly = TRUE, 	
		VacCampStrategies = "PROACTIVE", Efficacies_Start = 0.0, ImplementationDelays = 0, ImmunityDelays = 0, Efficacies_CamelControls = 0.5) 
Trees_CamelsOnly 	= read.table(file = file.path(CppOutputDirectory, paste0("Trees_CF_", ModelRun_CamelsOnly$OutputFolderNames, ".txt")), header = TRUE)

## Change H2Htype within Tree Data. DifferentRegion = 2 unhelpfully combines animal reservoir
## transmission and human-to-human transmission between regions.

AnimalReservoirLabelNum = 3
ChangeH2HLabel_AnimalTransmission = function(TreeData)
{
	TreeData$H2Htype[which(TreeData$H2Htype == DifferentRegion & TreeData$Infector == AnimalReservoir)] = AnimalReservoirLabelNum
	return(TreeData)
}

Trees_noCFs					= ChangeH2HLabel_AnimalTransmission(Trees_noCFs					)
Trees_Pro_ModEff_ShortLag	= ChangeH2HLabel_AnimalTransmission(Trees_Pro_ModEff_ShortLag	)
Trees_Re_ModEff_SlowReact	= ChangeH2HLabel_AnimalTransmission(Trees_Re_ModEff_SlowReact	)
Trees_CamelsOnly			= ChangeH2HLabel_AnimalTransmission(Trees_CamelsOnly			)

### want functions that: 
## i) 	extracts single transmission tree from larger data frame 
## ii) 	creates nodes, edges and graphs from single transmission tree. 
## iii) makes plots given nodes, edges etc. 

## i) extracts single transmission tree from larger data frame 
GetSingleTreeDataFrame = function(TreeData, iteration = NULL, CFnum = 0)
{
	if (is.null(iteration))	id.to.keep <- tail(TreeData$iteration, 1)	else  id.to.keep <- iteration
	TreeData = TreeData[TreeData$iteration == id.to.keep, ]
	if (!is.null(TreeData$CFnum))	#### i.e. if you are doing a CF tree
	{
		if ( (CFnum > max(TreeData$CFnum) | (CFnum < 0))) stop("GetSingleTreeDataFrame error: CFnum incorrect")
		if (is.null(CFnum)) CF_to_keep = tail(TreeData$CFnum, 1)  	else  CF_to_keep <- CFnum #### if you've not provided a particular CFnum, pick the last one as default. 
		TreeData = TreeData[TreeData$CFnum == CF_to_keep,]
	}
	return(TreeData)
}

## ii) creates nodes, edges and graphs from single transmission tree. 
Create_NodesEdgesGraph = function(TreeData = NULL)
{
	## We create an igraph using the edge list from the simulation.
	## Attributes for vertices will be added later on, after ensuring re-ordering of the nodes.
	
	Nodes 			<- with(TreeData, data.frame(id=CaseNumber, hospital=hospital))
	rownames(Nodes) <- Nodes$id
	Edges 			<- with(TreeData, data.frame(from = Infector, to = CaseNumber, type = H2Htype))
	
	## graph
	Graph	<- graph_from_data_frame(Edges, directed = TRUE)
	Nodes 	<- Nodes[get.vertex.attribute(Graph, "name"), ]	
	
	list(Graph = Graph, Edges = Edges, Nodes = Nodes)
}
## iii) makes plots given nodes, edges etc. 
GraphLayout = function(Graph, Method, graphIter)
{
	if (Method == "fr" ) Layout <- layout_with_fr			(Graph, niter	= graphIter)
	if (Method == "kk" ) Layout <- layout_with_kk			(Graph, maxiter	= graphIter)
	if (Method == "drl") Layout <- layout_with_drl			(Graph)
	if (Method == "rt" ) Layout <- layout.reingold.tilford	(Graph)
	
	return(Layout)
}

EdgeCols = c("orange", "green", "blue", "red")

Method = "rt"

MakePlotGivenGraph = function(GraphList, Method = "fr", graphIter = 5e4, IncTransmissionLegend = FALSE, 
		IncRegionLegend = FALSE, Layout = NULL, VertexLabel = "", FigLabel = "A", Title = "Title", ...)
{
	if (!(Method %in% c(c("fr", "kk", "drl", "rt")))) stop("MakePlotGivenGraph error: Method not recognized")
	### extract relevant parts of GraphList
	Graph <- GraphList$Graph
	Nodes <- GraphList$Nodes
	Edges <- GraphList$Edges
	## CHANGE ATTRIBUTES
	e.col 		<- EdgeCols[Edges$type + 1]
	col4edges 	<- unique(e.col)
	#print(col4edges)
	
	## vertex size = mean individual R
	Nodes$size 						<- 1.5
	Nodes$size[is.na(Nodes$size)] 	<- 1.5
	
	## vertex color = region
	Nodes$region <- 1
	
	Nodes$col <- factor(Nodes$region)
	K <- length(unique(Nodes$region))
	set.seed(1)
	reg.col <- sample(funky(K))
	levels(Nodes$col) <- reg.col
	Nodes$col <- as.character(Nodes$col)
	Nodes$col[is.na(Nodes$col)] <- transp("black")
	
	reg.annot <- levels(factor(Nodes$region))
	reg.annot[reg.annot == "Al Hudud ash Shamaliyah"] <- "Al Hudud ash S."
	
	set.seed(12)
	
	if (is.null(Layout)) Layout = GraphLayout(Graph, Method, graphIter)
	
	if (VertexLabel == "CaseNumber") VertexLabel = Nodes$id
	
	plot(Graph, vertex.color = Nodes$col, vertex.label = VertexLabel, 
			vertex.size = Nodes$size, edge.color = e.col, edge.width = 1, edge.arrow.mode = 0, layout = Layout, ...)
	
	if (IncTransmissionLegend	) legend("bottomright"	, 
				lwd = 3, lty = 1, col = EdgeCols, 
				legend = c("Within hospital", "Within region", "Between regions", "Animal reservoir"), 
				title = "Transmission type", cex = 1.3)
	if (IncRegionLegend			) legend("bottom"	, fill = reg.col, legend = reg.annot, title = "Region", ncol = 4, inset = c(0,-.12), cex = 1.5)
	
	if (!is.null(FigLabel))
	{
		par (xpd = TRUE)
		text(x = -1.25, y = 1.4, label = FigLabel, cex = 4)
	}
	text(x = 0, y = 1.35, label = Title, cex = 1.5, font = 2)
}
MakeSeveralPlotsGivenGraph 	= function(GraphList, Methods = c("fr", "kk", "drl", "rt")	, graphIter = 5e4, IncLegend = FALSE, Layout = NULL, VertexLabel = "")
{
	### set up par
	if (length(Methods) == 1) par(mfrow = c(2,2)) else
	if (length(Methods) == 2) par(mfrow = c(2,1)) else par(mfrow = c(2,2))
	
	for (Method in Methods)
		MakePlotGivenGraph(GraphList = GraphList, Method = Method, IncLegend = IncLegend)
}

# Tree 1
LastSim_No_CF_1 	<- GetSingleTreeDataFrame(Trees_noCFs)
G_List_No_CF_1 		<- Create_NodesEdgesGraph(TreeData = LastSim_No_CF_1)

LastSim_CF_1_3 		<- GetSingleTreeDataFrame(Trees_Re_ModEff_SlowReact)
G_List_CF_1_3 		<- Create_NodesEdgesGraph(TreeData = LastSim_CF_1_3)

LastSim_CF_2_1		<- GetSingleTreeDataFrame(Trees_Pro_ModEff_ShortLag, iteration = 5000)
G_List_CF_2_1 		<- Create_NodesEdgesGraph(TreeData = LastSim_CF_2_1)

LastSim_CF_2_4 		<- GetSingleTreeDataFrame(Trees_CamelsOnly, iteration = 3000)
G_List_CF_2_4 		<- Create_NodesEdgesGraph(TreeData = LastSim_CF_2_4)  # switched iteration here as legend easier to see on final one. 


png(file.path(PlotsDirectory, "IncPlotsAndTrees.png"), units = "in", res = 700, width = 15*2.5/3, height = 7.5)
#pdf(file.path(PlotsDirectory, "IncPlotsAndTrees.pdf"), width = 13*2.5/3, height = 2.5)
par(mar = c(1, 1, 1, 1))
LayoutMatrix = c(1,1,2,2)
LayoutMatrix = rbind(LayoutMatrix, c(3,4,5,6))
layout(LayoutMatrix, heights = c(1.5, 1))
par(xpd = TRUE)
par(mar= c(9.1, 7.1, 2.1, 0.8)) # no margins order: bottom, left, top, and right. The default is c(5.1, 4.1, 4.1, 2.1).

CasesCol = "grey70"              
DeathsCol = "firebrick1"             
DeathsCol = "black"             
HCW_Col = "blue"               
nHCW_Col = "green2"            
Deaths_HCW_Col = "gold"      
Deaths_nHCW_Col = "orangered"    

CEXLEG = 2

b <- barplot(AllCases_Inc, col = CasesCol, border = NA, width = 1, cex.axis = 1.7)
barplot(Cases_HCW, col = HCW_Col, border = NA, add = TRUE, width = 1, cex.axis = 1.7)
text(x = b[DatesForAxis_Int,1], y = -12, label = DateLabels, srt = 35, cex = 1.8) # have to place text at coords where barplot is actually drawn - because for some reason bar widths aren't equal to 1.
axis(1, at=b[DatesForAxis_Int,1], labels = FALSE)
legend("topleft", legend = c("all cases (n = 681)", "HCW cases (n = 187)"), col = c(CasesCol, HCW_Col), pch = 15, cex = CEXLEG)
text(x = -17.2, y = max(AllCases_Inc, na.rm = T) * 1.0, label = "A", cex = 4)

b <- barplot(AllDeaths_Inc, col = DeathsCol, border = NA, width = 1, cex.axis = 1.7)
barplot(Deaths_HCW, col = Deaths_HCW_Col, border = NA, add = TRUE, width = 1, cex.axis = 1.7)
text(x = b[DatesForAxis_Int,1], y = -4, label = DateLabels, srt = 35, cex = 1.8) # have to place text at coords where barplot is actually drawn - because for some reason bar widths aren't equal to 1.
axis(1, at=b[DatesForAxis_Int,1], labels = FALSE)
legend("topleft", legend = c("all deaths (n = 276)", "HCW deaths (n = 15)"), col = c(DeathsCol, Deaths_HCW_Col), pch = 15, cex = CEXLEG)
text(x = -17.2, y = max(AllDeaths_Inc, na.rm = T)* 1.0, label = "B", cex = 4)

par(mar = c(0, 0, 4.1, 0), xpd = TRUE)
CEXMAIN = 1.5
MakePlotGivenGraph(GraphList = G_List_No_CF_1, FigLabel = "C", Title = paste0("No vaccine:\n", dim(LastSim_No_CF_1)[1], " cases\n"), cex.main = 100,  Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = "", cex.main = CEXMAIN)
MakePlotGivenGraph(GraphList = G_List_CF_1_3 , FigLabel = "D", Title = paste0("Reactive (regional), VE = 60%,\n20-year mean duration,\n14-day react-time: ", dim(LastSim_CF_1_3)[1], " cases")	, Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = "", cex.main = CEXMAIN) 
MakePlotGivenGraph(GraphList = G_List_CF_2_1 , FigLabel = "E", Title = paste0("Proactive, VE = 60%,\n20-year mean duration,\n6-month wait: ", dim(LastSim_CF_2_1)[1], " cases")	, Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = "", cex.main = CEXMAIN)
MakePlotGivenGraph(GraphList = G_List_CF_2_4 , FigLabel = "F", Title = paste0("Camel control measures only,\n50% effective: ", dim(LastSim_CF_2_4)[1], " cases\n")		, Method = Method, IncTransmissionLegend = TRUE, IncRegionLegend = FALSE, VertexLabel = "", cex.main = CEXMAIN)
dev.off()




