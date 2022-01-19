
#### amend this script to use functions. Would be nice to say loop over a single variable (e.g. Active/Passive), without having to do all combinations of everything, For now though just 

rm(list=ls(all=TRUE)) 
options(width = 172L)
require(here)


ls()
ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppRootDirectory 	= file.path(ProjectDirectory, "MERS_Vac"	, "MERS_Vac") 
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
RawDataDirectory 	= file.path(ProjectDirectory, "Data"		)
PlotsDirectory		= file.path(ProjectDirectory, "Plots"		)
LocalParamFileDir	= file.path(ProjectDirectory, "ParamFiles"	)

#outputdirectory 	= LocalParamFileDir
outputdirectory 	= "\\\\fi--didenas1\\dengue\\Danny\\MERS\\"
#outputdirectory 	= "\\\\fi--didenas1\\dengue\\Danny\\MERS"

WriteParamFiles 		= TRUE
WriteSingleBatchFiles 	= TRUE
WriteManyBatchFiles	  	= TRUE
ParamFileNames 			= ""

source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))

### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 
### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 
### === 		Make "Single" Batch file. 
### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 

#"job submit /scheduler:fi--didemrchnb /jobtemplate:32Core 		/numnodes:1 /singlenode:false /workdir:\\fi--san02\homes\cpm14\Tampines6 VectorModel6.exe %1 %2 %3 %4"
#"job submit /scheduler:fi--didemrchnb /jobtemplate:GeneralNodes /numnodes:1 /singlenode:false /workdir:\\fi--didenas1-app\dengue\Danny /stdout:stdout_AS_PRIME_%1 /stderr:stderr_AS_PRIME_%1 D_MCMC.exe %1"

WriteSingleBatchFiles = function(ManyANDSingleFileName, Nodes = "GeneralNodes", stderrlabel = "_%1", stdoutlabel = "_%1", 
		Whole_Node = TRUE, NumCores = NULL, Exe_name = "MERS_Vac.exe", Exe_args = "%1", BatchOriginDir = "/workdir:\\\\fi--didenas1-app\\dengue\\Danny\\MERS", WRITE_FILES = TRUE, 
		AddBatchNameToOutputStreams = TRUE)
{
	if (Exe_name != "MERS_Vac.exe") warning("Exe_name != MERS_Vac.exe")
	if (ManyANDSingleFileName == "") stop("Must Provide ManyANDSingleFileName argument")
	
	SingleBatchFileName = paste0(ManyANDSingleFileName, "_Single.bat")
	
	if (Whole_Node) NumNodes_SingleNode_string = " /numnodes:1 /singlenode:false" 		else
					NumNodes_SingleNode_string = paste0(" /numcores:", NumCores, "-", NumCores, " /singlenode:true") 			#### check this with Wes - you don't actually know the pattern you're just cargo culting it. 
	ClusterString 	= paste0("job submit /scheduler:fi--didemrchnb /jobtemplate:", Nodes, NumNodes_SingleNode_string)
	
	if (AddBatchNameToOutputStreams) stderr_BatchAddOn = ManyANDSingleFileName else stderr_BatchAddOn = ""
	if (AddBatchNameToOutputStreams) stdout_BatchAddOn = ManyANDSingleFileName else stdout_BatchAddOn = ""
	
	stderrlabel = paste0("/stderr:STDERROUTs\\stderr_", stderr_BatchAddOn, stderrlabel) ### add output stuff Wes told you to add. 
	stdoutlabel = paste0("/stdout:STDERROUTs\\stdout_", stdout_BatchAddOn, stdoutlabel)
	
	LongString_single = paste(ClusterString, BatchOriginDir, stdoutlabel, stderrlabel, Exe_name, Exe_args)
#	print(LongString_single)
	if (WRITE_FILES) write.table(LongString_single, file = file.path(outputdirectory, SingleBatchFileName), row.names = F, col.names = F, quote = F, sep = "\t")
}

### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 
### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 
### === 		Make Param Files and "Many" Batch files that launches multiple jobs 
### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 



options(scipen=999)

WriteParamFile = function(ModelRun, NumIterations = "11000", BurnIn = "1000", WRITE_FILES = TRUE, WarningMessage = NULL)
{
	#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
	#### #### #### #### #### Error-trapping / Checks
	#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
	
	if (!is.null(WarningMessage)) warning(WarningMessage)
	
	if (class(BurnIn) 						!= "character") stop("BurnIn must be character"						)
	if (class(NumIterations) 				!= "character") stop("NumIterations must be a character"				)
#	if (class(seed1						) 	!= "character") stop("seed1 must be character"						)
#	if (class(seed2						) 	!= "character") stop("seed2 must be character"						)
	if (as.numeric(BurnIn) >= as.numeric( NumIterations))	stop("BurnIn >= NumIterations"						)
	
	ParamFileName = ChooseOutputString(ModelRun, StringPrefix = "Params", Folder = FALSE) 
	
	### add .txt extension
	ParamFileName 	= paste0(ParamFileName, ".txt")
	ParamTable 		= rbind(t(ModelRun), NumIterations, BurnIn)
	ParamFileName 	<<- ParamFileName

	class(ParamTable)
	if (!(any(ParamFileName == ParamFileNames)))
	{
		ParamTable <<- ParamTable
		
		if (WRITE_FILES) write.table(ParamTable, file = file.path(outputdirectory, ParamFileName),	row.names = T, col.names = F, quote = F, sep = "\t")
		
		#### add to ParamFileNames (globally)
		ParamFileNames <<- c(ParamFileNames, ParamFileName )
		
		return(ParamFileName)
	}
}

WriteMultiBatchFile = function(ParamFileNames, ManyANDSingleFileName = ManySingleFileName)
{
	if (ParamFileNames[1] == "") ParamFileNames = ParamFileNames[-1]
	ManyBatchFileName 		= paste0(ManyANDSingleFileName, "_Many.bat")
	SingleBatchFileName 	= paste0(ManyANDSingleFileName, "_Single.bat")
	
	ManyBatchFileName <<- ManyBatchFileName
	SingleBatchFileName <<- SingleBatchFileName
	
	MultiJobString = ""
	for (ParamFileName in ParamFileNames)
	{
		CallString = paste("call", SingleBatchFileName, ParamFileName, "\n")
		MultiJobString = paste0(MultiJobString, CallString)
	}
	#MultiJobString <<- MultiJobString
	write.table(MultiJobString, file = file.path(outputdirectory, ManyBatchFileName),	row.names = F, col.names = F, quote = F, sep = "\t")
}

#source(paste0(R_ScriptDirectory, "DirectoryFunctions2.R"))

#ManySingleFileName = "FirstClusterRuns"	### 
#ManySingleFileName = "SecondClusterRuns"	### 
#ManySingleFileName = "ThirdClusterRuns"	### 
#ManySingleFileName = "FourthClusterRuns"	### changed number of times trees are outputted.  
#ManySingleFileName = "FifthClusterRuns"	### Additional efficacy values considered  
#ManySingleFileName = "FifthClusterRuns_fillinblanks"	### Additional efficacy values considered  

### reactive
#ManySingleFileName = "Reactive"	
#ManySingleFileName = "Reactive_wRegionalAndCamels"	
##ModelRuns = DefineModelRuns(ImplementationDelays = seq(0,28,4), ImmunityDelays = 14, Efficacies_Start = seq(0.05, 1, by = 0.05))
#ModelRuns = DefineModelRuns(#NewRunsOnly	= TRUE, 
#		ImplementationDelays = seq(0,28,2), ImmunityDelays = 14, Efficacies_Start = seq(0.05, 1, by = 0.05), 
#		ReactiveAtHospitalLevel_And_Not = 1:0, Vaccinate_Camels_And_Not = 0:1)

#ManySingleFileName = "Reactive_RegionalTest"
#ModelRuns = DefineModelRuns(
#		#NewRunsOnly	= TRUE, 
#		ImplementationDelays = seq(0,28,14), ImmunityDelays = 14, Efficacies_Start = seq(0.05, 1, by = 0.05), 
#		ReactiveAtHospitalLevel_And_Not	= 1, 
#		Vaccinate_Camels_And_Not 		= 0)

#ManySingleFileName = "Reactive_CamelsTest"
#ModelRuns = DefineModelRuns(
#		#NewRunsOnly	= TRUE, 
#		ImplementationDelays = seq(0,28,14), ImmunityDelays = 14, Efficacies_Start = seq(0.05, 1, by = 0.05), 
#		ReactiveAtHospitalLevel_And_Not	= 0, 
#		Vaccinate_Camels_And_Not 		= 1)


### proactive
#ManySingleFileName = "Proactive"	### Additional efficacy values considered  
#ModelRuns = DefineModelRuns(#NewRunsOnly = TRUE,
#		VacCampStrategies = "PROACTIVE", ImplementationDelays = 0, ImmunityDelays = 0, ##  ImplementationDelays = 0 and ImmunityDelays = 0 are an approximation.
#		Efficacies_Start = seq(0.05, 1, by = 0.05), VaccineDurations = c(0, 20, 10, 5, 2, 1), TimesSinceVaccination = c(0.5, 1:10))

#ManySingleFileName = "ProactiveCamels"	### Additional efficacy values considered  
#ModelRuns = DefineModelRuns(#NewRunsOnly = TRUE,
#		Vaccinate_Camels_And_Not = 1,
#		VacCampStrategies = "PROACTIVE", ImplementationDelays = 0, ImmunityDelays = 0, ##  ImplementationDelays = 0 and ImmunityDelays = 0 are an approximation.
#		Efficacies_Start = seq(0.05, 1, by = 0.05), VaccineDurations = c(0, 20, 10, 5, 2, 1), TimesSinceVaccination = c(0.5, 1:10))

#ManySingleFileName = "ReactiveHosp_FillInDeathBlanks"	
#ModelRuns = DefineModelRuns(ImplementationDelays = seq(0,28,4), ImmunityDelays = 14, Efficacies_Start = seq(0.05, 1, by = 0.05))
#ModelRuns = DefineModelRuns(#NewRunsOnly	= TRUE, 
#		ImplementationDelays = seq(0,28,2), ImmunityDelays = 14, Efficacies_Start = seq(0.05, 1, by = 0.05), 
#		ReactiveAtHospitalLevel_And_Not = 1, Vaccinate_Camels_And_Not = 0)

## reactive at national level
#ManySingleFileName = "Reactive_wNationalAndCamels"	
##ModelRuns = DefineModelRuns(ImplementationDelays = seq(0,28,4), ImmunityDelays = 14, Efficacies_Start = seq(0.05, 1, by = 0.05))
#ModelRuns = DefineModelRuns(#NewRunsOnly	= TRUE, 
#		ImplementationDelays = seq(0,28,2), ImmunityDelays = 14, Efficacies_Start = seq(0.05, 1, by = 0.05), 
#		ReactLevels = "NATIONAL", Vaccinate_Camels_And_Not = 0:1)

#ManySingleFileName = "ReactiveHigherCamelCoverage_60Percent"	
##ModelRuns = DefineModelRuns(ImplementationDelays = seq(0,28,4), ImmunityDelays = 14, Efficacies_Start = seq(0.05, 1, by = 0.05))
#ModelRuns = DefineModelRuns(#NewRunsOnly	= TRUE, 
#		ImplementationDelays = seq(0,28,2), ImmunityDelays = 14, Efficacies_Start = seq(0.05, 1, by = 0.05), 
#		ReactLevels = c("HOSPITAL", "REGIONAL", "NATIONAL"), Vaccinate_Camels_And_Not = 1, Coverages_Camels = 0.6)

#ManySingleFileName = "ProactiveHigherCamelCoverage_60Percent"
#ModelRuns = DefineModelRuns(NewRunsOnly = TRUE,
#		Vaccinate_Camels_And_Not = 1,
#		VacCampStrategies = "PROACTIVE", ImplementationDelays = 0, ImmunityDelays = 0, ##  ImplementationDelays = 0 and ImmunityDelays = 0 are an approximation.
#		Efficacies_Start = seq(0.05, 1, by = 0.05), VaccineDurations = c(0, 20, 10, 5, 2, 1), TimesSinceVaccination = c(0.5, 1:10), Coverages_Camels = 0.6)

## camel control measures
# proactive, or no humans, want "efficacies" of 0.1, 0.2, 0.3, 0.4, 0.5.

#ManySingleFileName = "Proactive_WithCamels"
#ModelRuns = DefineModelRuns(NewRunsOnly = TRUE,
#		Efficacies_CamelControls = seq(0, 0.5, 0.1), ## 0 here refers to no camels, or equivalently camel interventions having zero efficacy.
#		VacCampStrategies = "PROACTIVE", ImplementationDelays = 0, ImmunityDelays = 0, ##  ImplementationDelays = 0 and ImmunityDelays = 0 are an approximation.
#		Efficacies_Start = seq(0.00, 1, by = 0.05), ## 0 here refers to no humans, or equivalently human vaccines having zero efficacy.
#		VaccineDurations = c(0, 20, 15, 10, 5, 2, 1), TimesSinceVaccination = c(0.5, 1:10))

ManySingleFileName = "Reactive_WithCamels"	
ModelRuns = DefineModelRuns(
		NewRunsOnly = FALSE, 
		VacCampStrategies = "REACTIVE", 
		Efficacies_CamelControls = seq(0.1, 0.5, 0.1),
		ImplementationDelays = seq(0,28,2), 
		ImmunityDelays = 14, 
		Efficacies_Start = seq(0.05, 1, by = 0.05), 
		ReactLevels	= c("HOSPITAL", "REGIONAL", "NATIONAL")							, #c("HOSPITAL", "REGIONAL", "NATIONAL")
		)

#ManySingleFileName = "MissedRuns_211110"	
#ModelRuns = read.table(file = file.path(CppOutputDirectory, paste0("MissedRuns_211110.txt")), header = T, sep = "\t")
		

dim(ModelRuns)
str(ModelRuns)

#dim(ModelRunsBeforeCull)
#str(ModelRunsBeforeCull)



#VaccineDurations						= 20000.0						,
#TimesSinceVaccination					= 0								,


#ModelRuns$OutputFolderNames
length(ModelRuns$OutputFolderNames) == length(unique(ModelRuns$OutputFolderNames))
ModelRun = ModelRuns[1,]

WriteSingleBatchFiles(ManyANDSingleFileName = ManySingleFileName, Whole_Node = FALSE, NumCores = 1)

MR_index = 1
for (MR_index in 1:dim(ModelRuns)[1]) WriteParamFile(ModelRun = ModelRuns[MR_index,])

WriteMultiBatchFile(ParamFileNames)

warnings()








