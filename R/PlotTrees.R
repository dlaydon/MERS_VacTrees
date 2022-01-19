#### This script contains your changes to script plotATree.1.5.R, (the original version of which is in Dropbox (SPH Imperial College)\nCoV\KSA\simon_EM\transmissionTree\hugATree, but which you've already made some small changes to)
#### You want to turn code into function that will enable you to plot many trees

rm(list = ls())
ls()

options(width = 141L)

## LOAD PACKAGES ##
#install.packages(c("igraph", "adegenet", "RColorBrewer"))
library(igraph)
library(adegenet)
library(RColorBrewer)
require(here)
library(grid)
library(gridExtra)

## Want 2x5 figure showing multiple who-infected-whom trees.
# 1st column: no vaccine; 
# 2nd column: proactive moderate efficacy & duration, short lag; 
# 3rd column: proactive moderate efficacy & duration, long lag; 
# 4th column: reactive with e.g. 10-day delay at regional level. 
# 5th column: Camels vaccinated (say Eff = 30%)
# 2nd row: alternative tree to illustrate that they are modelled, not measured.

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

here("R")
here("MERS_Vac/MERS_Vac")

source(here("R/DirectoryFunctions.R"))
source(here("R/LayoutPlottingFunctions.R")) 

# 1st column: no vaccine; 
Trees_noCFs 	<- read.table(here("Output/ExampleTrees_nonCF.txt"), header = TRUE)
str(Trees_noCFs)
head(Trees_noCFs)

# 2nd column: proactive moderate efficacy & duration, short lag; 
ModelRun_Pro_ModEff_ShortLag 	= DefineModelRuns(IncludeCompletedRunsOnly = TRUE, VacCampStrategies = "PROACTIVE", Efficacies_Start = 0.6, VaccineDurations = 5, TimesSinceVaccination = 0.5) 
Trees_Pro_ModEff_ShortLag 		= read.table(file = file.path(CppOutputDirectory, paste0("Trees_CF_", ModelRun_Pro_ModEff_ShortLag$OutputFolderNames, ".txt")), header = TRUE)

# 3rd column: proactive moderate efficacy & duration, long lag; 
ModelRun_Pro_ModEff_LongLag = DefineModelRuns(IncludeCompletedRunsOnly = TRUE, 	
		VacCampStrategies = "PROACTIVE", Efficacies_Start = 0.6, VaccineDurations = 20, TimesSinceVaccination = 5) 
Trees_Pro_ModEff_LongLag 	= read.table(file = file.path(CppOutputDirectory, paste0("Trees_CF_", ModelRun_Pro_ModEff_LongLag$OutputFolderNames, ".txt")), header = TRUE)

# 4th column: reactive with e.g. 10-day delay at regional level. 
ModelRun_Re_ModEff_SlowReact = DefineModelRuns(IncludeCompletedRunsOnly = TRUE, 	
		VacCampStrategies = "REACTIVE", ReactLevels = "REGIONAL", Efficacies_Start = 0.6, ImplementationDelays = 14, ImmunityDelays = 14) 
Trees_Re_ModEff_SlowReact 	= read.table(file = file.path(CppOutputDirectory, paste0("Trees_CF_", ModelRun_Re_ModEff_SlowReact$OutputFolderNames, ".txt")), header = TRUE)

# 5th column: reactive with e.g. 10-day delay at regional level. 
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
Trees_Pro_ModEff_LongLag	= ChangeH2HLabel_AnimalTransmission(Trees_Pro_ModEff_LongLag	)
Trees_Re_ModEff_SlowReact	= ChangeH2HLabel_AnimalTransmission(Trees_Re_ModEff_SlowReact	)
Trees_CamelsOnly			= ChangeH2HLabel_AnimalTransmission(Trees_CamelsOnly			)
#range(Trees_CamelsOnly$iteration)

signif(table(Trees_noCFs$H2Htype) / dim(Trees_noCFs)[1], 3)
signif(table(Trees_CamelsOnly$H2Htype) / dim(Trees_CamelsOnly)[1], 3)
signif(table(Trees_Pro_ModEff_ShortLag$H2Htype) / dim(Trees_Pro_ModEff_ShortLag)[1], 3)

sum(table(Trees_noCFs$H2Htype) / dim(Trees_noCFs)[1])
sum(table(Trees_Pro_ModEff_ShortLag$H2Htype) / dim(Trees_Pro_ModEff_ShortLag)[1])

## HOSPITAL INFO
hospreg <- read.csv(here("Data/hospital-region.csv"))
head(hospreg)
regnames <- read.table(here("Data/regionID-names.txt"), row.names = 1)
regnames
regions <- regnames[as.character(hospreg$region.ID),]
names(regions) <- hospreg$hospital.ID
str(Trees_noCFs)

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

EdgeCols = c("yellow", "red", "blue", "green")
EdgeCols = c("orange", "green", "blue", "red")
#for (i in 1:dim(Edges)[1]) e.col[i] = EdgeCols[Edges$type[i]]
#
#e.col = EdgeCols[Edges$type]

GraphList = G_List_No_CF_1
MakePlotGivenGraph = function(GraphList, Method = "fr", graphIter = 5e4, IncTransmissionLegend = FALSE, 
		IncRegionLegend = FALSE, Layout = NULL, VertexLabel = "", ...)
{
	if (!(Method %in% c(c("fr", "kk", "drl", "rt")))) stop("MakePlotGivenGraph error: Method not recognized")
	### extract relevant parts of GraphList
	Graph <- GraphList$Graph
	Nodes <- GraphList$Nodes
	Edges <- GraphList$Edges
	## CHANGE ATTRIBUTES
#	e.col 		<- transp(fac2col(Edges$type + 1, col.pal = seasun))
	e.col 		<- EdgeCols[Edges$type + 1]
	col4edges 	<- unique(e.col)
	#print(col4edges)
	
	## vertex size = mean individual R
	Nodes$size 						<- 1.5 * (1 + (table(Edges$from)[rownames(Nodes)])^0.5)
	Nodes$size[is.na(Nodes$size)] 	<- 1.5
	
	## vertex color = region
	Nodes$region <- regions[as.character(Nodes$hospital)]
	
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
	Layout <<- Layout
	
	if (VertexLabel == "CaseNumber") VertexLabel = Nodes$id
	
	###  better to set par outside function
	#par(mar=c(4.5,.1,.1,.1))
	#par(mar=rep(0,4))
	#if (IncTransmissionLegend & IncRegionLegend) par(mar=c(4.5,.1,.1,.1)) else par(mar=rep(0,4))
	plot(Graph, vertex.color = Nodes$col, vertex.label = VertexLabel, vertex.size = Nodes$size, edge.color = e.col, edge.width = 1, edge.arrow.mode = 0, layout = Layout, ...)
	
	if (IncTransmissionLegend	) legend("bottomright"	, 
				lwd = 3, lty = 1, col = EdgeCols, legend = c("Within hospital", "Within region", "Between regions", "Animal reservoir"), 
				title = "Transmission type", cex = 1.15)
	if (IncRegionLegend			) legend("bottom"	, fill = reg.col, legend = reg.annot, title = "Region", ncol = 4, inset = c(0,-.12), cex = 1.5)
}
MakeSeveralPlotsGivenGraph 	= function(GraphList, Methods = c("fr", "kk", "drl", "rt")	, graphIter = 5e4, IncLegend = FALSE, Layout = NULL, VertexLabel = "")
{
	#### Want this function to display either: i) single tree with several methods; or ii) several trees with the same method. Right now it does i) only. 
	### set up par
	if (length(Methods) == 1) par(mfrow = c(2,2)) else
	if (length(Methods) == 2) par(mfrow = c(2,1)) else par(mfrow = c(2,2))
	
	for (Method in Methods)
		MakePlotGivenGraph(GraphList = GraphList, Method = Method, IncLegend = IncLegend)
}

# ROW 1
LastSim_No_CF_1 	<- GetSingleTreeDataFrame(Trees_noCFs)
G_List_No_CF_1 		<- Create_NodesEdgesGraph(TreeData = LastSim_No_CF_1)

LastSim_CF_1_1		<- GetSingleTreeDataFrame(Trees_Pro_ModEff_ShortLag)
G_List_CF_1_1 		<- Create_NodesEdgesGraph(TreeData = LastSim_CF_1_1)

LastSim_CF_1_2 		<- GetSingleTreeDataFrame(Trees_Pro_ModEff_LongLag)
G_List_CF_1_2 		<- Create_NodesEdgesGraph(TreeData = LastSim_CF_1_2)

LastSim_CF_1_3 		<- GetSingleTreeDataFrame(Trees_Re_ModEff_SlowReact)
G_List_CF_1_3 		<- Create_NodesEdgesGraph(TreeData = LastSim_CF_1_3)

LastSim_CF_1_4 		<- GetSingleTreeDataFrame(Trees_CamelsOnly, iteration = 5000) # switched iteration here as legend easier to see on final one. 
G_List_CF_1_4 		<- Create_NodesEdgesGraph(TreeData = LastSim_CF_1_4)

# ROW 2
LastSim_No_CF_2 	<- GetSingleTreeDataFrame(Trees_noCFs, iteration = 5000)
G_List_No_CF_2 		<- Create_NodesEdgesGraph(TreeData = LastSim_No_CF_2)

LastSim_CF_2_1		<- GetSingleTreeDataFrame(Trees_Pro_ModEff_ShortLag, iteration = 5000)
G_List_CF_2_1 		<- Create_NodesEdgesGraph(TreeData = LastSim_CF_2_1)

LastSim_CF_2_2 		<- GetSingleTreeDataFrame(Trees_Pro_ModEff_LongLag, iteration = 5000)
G_List_CF_2_2 		<- Create_NodesEdgesGraph(TreeData = LastSim_CF_2_2)

LastSim_CF_2_3 		<- GetSingleTreeDataFrame(Trees_Re_ModEff_SlowReact, iteration = 5000)
G_List_CF_2_3 		<- Create_NodesEdgesGraph(TreeData = LastSim_CF_2_3)

LastSim_CF_2_4 		<- GetSingleTreeDataFrame(Trees_CamelsOnly)
G_List_CF_2_4 		<- Create_NodesEdgesGraph(TreeData = LastSim_CF_2_4)  # switched iteration here as legend easier to see on final one. 


Method = "fr"	### takes forever
Method = "rt"
#Method = "kk"
#layout(t(c(1,2,3,4)))
#par(mar=rep(0,4))
#par(mar=rep(0,4), mfrow = c(2,4))
png(file.path(PlotsDirectory, "ExampleTrees.png"), units = "in", res = 700, width = 10*2.5/3, height = 2.5)
par(mfrow = c(1,4), mar = c(0, 0, 3.1, 0))

## Want 2x4 figure showing multiple who-infected-whom trees.
# 1st column: no vaccine; 
# 2nd column: proactive moderate efficacy & duration, short lag; 
# 3rd column: proactive moderate efficacy & duration, long lag; 
# 4th column: reactive with e.g. 10-day delay at regional level. 
# 2nd row: alternative tree to illustrate that they are modelled, not measured.

MakePlotGivenGraph(GraphList = G_List_No_CF_1, main = paste0("(A) No vaccine,\ntree A: ", dim(LastSim_No_CF_1)[1], " cases\n"), cex.main = 100,  Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = "")
#MakePlotGivenGraph(GraphList = G_List_CF_1_1 , main = paste0("PROACTIVE, 60% efficacy,\n5-year mean duration, 6-month wait,\nsample A: ", dim(LastSim_CF_1_1)[1], " cases")	, Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = ""		)
#MakePlotGivenGraph(GraphList = G_List_CF_1_2 , main = paste0("PROACTIVE, 60% efficacy,\n5-year mean duration, 5-year lag,\nsample A: ", dim(LastSim_CF_1_2)[1], " cases")	, Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = ""		)
#MakePlotGivenGraph(GraphList = G_List_CF_1_3 , main = paste0("REACTIVE to regional outbreak,\n60% efficacy, 10-day delay,\nsample A: ", dim(LastSim_CF_1_3)[1], " cases")	, Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = ""	) 
#MakePlotGivenGraph(GraphList = G_List_CF_1_4 , main = paste0("ONLY CAMELS considered,\ncontrol measures 50% effective,\nsample A: ", dim(LastSim_CF_1_4)[1], " cases")		, Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = ""	) 

MakePlotGivenGraph(GraphList = G_List_No_CF_2, main = paste0("(B) No vaccine,\ntree B: ", dim(LastSim_No_CF_2)[1], " cases\n")													, Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = "")
MakePlotGivenGraph(GraphList = G_List_CF_2_1 , main = paste0("(C) Proactive, 60% efficacy,\n10-year mean duration,\n6-month wait, tree A: ", dim(LastSim_CF_2_1)[1], " cases")	, Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = ""		)
#MakePlotGivenGraph(GraphList = G_List_CF_2_2 , main = paste0("PROACTIVE, 60% efficacy,\n5-year mean duration, 5-year lag,\nsample B: ", dim(LastSim_CF_2_2)[1], " cases")	, Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = ""		)
#MakePlotGivenGraph(GraphList = G_List_CF_2_3 , main = paste0("REACTIVE to regional outbreak,\n60% efficacy, 10-day delay,\nsample B: ", dim(LastSim_CF_2_3)[1], " cases")	, Method = Method, IncTransmissionLegend = FALSE, IncRegionLegend = FALSE, VertexLabel = ""	) 
MakePlotGivenGraph(GraphList = G_List_CF_2_4 , main = paste0("(D) Camel control measures\nonly, 50% effective,\ntree A: ", dim(LastSim_CF_2_4)[1], " cases")		, Method = Method, IncTransmissionLegend = TRUE, IncRegionLegend = FALSE, VertexLabel = ""	) 

dev.off()



 

#MakePlotGivenGraph(GraphList = G_List_CF4	, Method = Method	, IncLegend = FALSE, VertexLabel = "", main = paste0("50% coverage, 50% efficacy:\n"	, dim(LastSim_CF4)[1], " cases")) 
#MakePlotGivenGraph(GraphList = G_List_CF5	, Method = Method	, IncLegend = FALSE, VertexLabel = "", main = paste0("50% coverage, 50% efficacy:\nHealthcare workers only: "	, dim(LastSim_CF5)[1], " cases")) 

dev.off()


##### try to remove all hospital transmission here = it doesn't work (or at least the tree's look mental)
#LastSim_CF5 	<- LastSim[LastSim$H2Htype != 0, ]
#G_List_CF5	<- Create_NodesEdgesGraph(TreeData = LastSim_CF5)



#par(mar=rep(0,4), mfrow=c(2,2))				## bottom, left, top, and right
par(mar = c(0,0,2,0), mfrow = c(2,3))		## bottom, left, top, and right

#MakeSeveralPlotsGivenGraph	(GraphList = G_List)


#MakePlotGivenGraph			(GraphList = G_List, Method = "fr"	, IncLegend = FALSE)
#MakePlotGivenGraph			(GraphList = G_List, Method = "kk"	, IncLegend = FALSE)
#MakePlotGivenGraph			(GraphList = G_List, Method = "drl"	, IncLegend = FALSE)


#MakePlotGivenGraph			(GraphList = G_List		, Method = Method	, IncLegend = FALSE, VertexLabel = "CaseNumber")
#MakePlotGivenGraph			(GraphList = G_List_CF0	, Method = Method	, IncLegend = FALSE, VertexLabel = "CaseNumber")
#MakePlotGivenGraph			(GraphList = G_List_CF1	, Method = Method	, IncLegend = FALSE, VertexLabel = "CaseNumber")

#MakePlotGivenGraph			(GraphList = G_List		, Method = Method	, IncLegend = FALSE, VertexLabel = "", main = paste0("No vaccine:\n", dim(LastSim)[1], " cases"))
identical(G_List, G_List_CF2)

png(file = paste0(ProjectDirectory, "Talks", slash(), "SimonTree.png"), res = 300, height = 10, width = 10, units = "in")
#par(mar = c(0,0,1,0), mfrow=c(1,1))
par(mar = c(0,0,0,0), mfrow=c(1,1))
MakePlotGivenGraph			(GraphList = G_List		, Method = Method)
dev.off()

png(file = paste0(ProjectDirectory, "Talks", slash(), "VitalNodesPruned.png"), res = 300, height = 10, width = 10, units = "in")
#par(mar = c(0,0,1,0), mfrow=c(1,1))
par(mar = c(0,0,0,0), mfrow=c(1,1))
MakePlotGivenGraph			(GraphList = G_List_CF0	, Method = Method)
dev.off()

png(file = paste0(ProjectDirectory, "Talks", slash(), "FullCoverage_LoEff.png"), res = 300, height = 10, width = 10, units = "in")
#par(mar = c(0,0,1,0), mfrow=c(1,1))
par(mar = c(0,0,0,0), mfrow=c(1,1))
MakePlotGivenGraph			(GraphList = G_List_CF1	, Method = Method)
dev.off()

png(file = paste0(ProjectDirectory, "Talks", slash(), "FullCoverage_HiEff.png"), res = 300, height = 10, width = 10, units = "in")
#par(mar = c(0,0,1,0), mfrow=c(1,1))
par(mar = c(0,0,0,0), mfrow=c(1,1))
MakePlotGivenGraph			(GraphList = G_List_CF2	, Method = Method)
dev.off()

png(file = paste0(ProjectDirectory, "Talks", slash(), "HalfCoverage_LoEff.png"), res = 300, height = 10, width = 10, units = "in")
#par(mar = c(0,0,1,0), mfrow=c(1,1))
par(mar = c(0,0,0,0), mfrow=c(1,1))
MakePlotGivenGraph			(GraphList = G_List_CF3	, Method = Method)
dev.off()

png(file = paste0(ProjectDirectory, "Talks", slash(), "HalfCoverage_HiEff.png"), res = 300, height = 10, width = 10, units = "in")
#par(mar = c(0,0,1,0), mfrow=c(1,1))
par(mar = c(0,0,0,0), mfrow=c(1,1))
MakePlotGivenGraph			(GraphList = G_List_CF4	, Method = Method)
dev.off()

png(file = paste0(ProjectDirectory, "Talks", slash(), "HalfCoverage_HiEff_HCW_only.png"), res = 300, height = 10, width = 10, units = "in")
#par(mar = c(0,0,1,0), mfrow=c(1,1))
par(mar = c(0,0,0,0), mfrow=c(1,1))
MakePlotGivenGraph			(GraphList = G_List_CF5	, Method = Method)
dev.off()



#MakePlotGivenGraph			(GraphList = G_List		, Method = Method	, IncLegend = FALSE, VertexLabel = "", main = paste0("No vaccine:\n", dim(LastSim)[1], " cases"))
MakePlotGivenGraph			(GraphList = G_List_CF0	, Method = Method	, IncLegend = FALSE, VertexLabel = "", main = "")
MakePlotGivenGraph			(GraphList = G_List_CF1	, Method = Method	, IncLegend = FALSE, VertexLabel = "", main = paste0("10% efficacy, 100% coverage:\n", dim(LastSim_CF1)[1], " cases"))
MakePlotGivenGraph			(GraphList = G_List_CF2	, Method = Method	, IncLegend = FALSE, VertexLabel = "", main = paste0("50% efficacy, 100% coverage:\n", dim(LastSim_CF2)[1], " cases"))
MakePlotGivenGraph			(GraphList = G_List_CF3	, Method = Method	, IncLegend = FALSE, VertexLabel = "", main = paste0("50% efficacy, 50% coverage:\n"	, dim(LastSim_CF3)[1], " cases")) 
MakePlotGivenGraph			(GraphList = G_List_CF4	, Method = Method	, IncLegend = FALSE, VertexLabel = "", main = paste0("50% efficacy, 50% coverage:\n"	, dim(LastSim_CF4)[1], " cases")) 
MakePlotGivenGraph			(GraphList = G_List_CF5	, Method = Method	, IncLegend = FALSE, VertexLabel = "", main = paste0("50% efficacy, 50% coverage:\nHealthcare workers only: "	, dim(LastSim_CF5)[1], " cases")) 








































