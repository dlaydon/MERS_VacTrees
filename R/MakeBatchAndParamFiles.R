
require(here)


ls()
ProjectDirectory 	= here()
R_ScriptDirectory 	= file.path(ProjectDirectory, "R"			)
CppRootDirectory 	= file.path(ProjectDirectory, "MERS_Vac"	, "MERS_Vac") 
CppOutputDirectory 	= file.path(ProjectDirectory, "Output"		) 
RawDataDirectory 	= file.path(ProjectDirectory, "Data"		)
PlotsDirectory		= file.path(ProjectDirectory, "Plots"		)
LocalParamFileDir	= file.path(ProjectDirectory, "ParamFiles"	)

outputdirectory 	= LocalParamFileDir

WriteParamFiles 		= TRUE
WriteSingleBatchFiles 	= TRUE
WriteManyBatchFiles	  	= TRUE
ParamFileNames 			= ""

source(file.path(R_ScriptDirectory, "DirectoryFunctions.R"))

### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 
### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 
### === 		Make "Single" Batch file. 
### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === ### === 

WriteSingleBatchFiles = function(ManyANDSingleFileName, Exe_name = "MERS_Vac.exe", Exe_args = "%1", WRITE_FILES = TRUE)
{
	if (Exe_name != "MERS_Vac.exe") warning("Exe_name != MERS_Vac.exe")
	if (ManyANDSingleFileName == "") stop("Must Provide ManyANDSingleFileName argument")
	
	SingleBatchFileName = paste0(ManyANDSingleFileName, "_Single.bat")
	LongString_single 	= paste(Exe_name, Exe_args)
	if (WRITE_FILES) write.table(LongString_single, file = file.path(outputdirectory, SingleBatchFileName), row.names = F, col.names = F, quote = F, sep = "\t")
}

options(scipen=999)

WriteParamFile = function(ModelRun, NumIterations = "11000", BurnIn = "1000", WRITE_FILES = TRUE, WarningMessage = NULL)
{
	#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
	#### #### #### #### #### Error-trapping / Checks
	#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
	
	if (!is.null(WarningMessage)) warning(WarningMessage)
	
	if (class(BurnIn) 						!= "character") stop("BurnIn must be character"						)
	if (class(NumIterations) 				!= "character") stop("NumIterations must be a character"				)
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
	
	ManyBatchFileName 	<<- ManyBatchFileName
	SingleBatchFileName <<- SingleBatchFileName
	
	MultiJobString = ""
	for (ParamFileName in ParamFileNames)
	{
		CallString = paste(SingleBatchFileName, ParamFileName, "\n")
		MultiJobString = paste0(MultiJobString, CallString)
	}
	write.table(MultiJobString, file = file.path(outputdirectory, ManyBatchFileName),	row.names = F, col.names = F, quote = F, sep = "\t")
}


## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ 
## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ 
## ** -- ^^ MAKE BATCH FILES
## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ 
## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ 

## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ 
## ** -- ^^ PROACTIVE

ManySingleFileName = "Proactive_WithCamels"
ModelRuns = DefineModelRuns(
		Efficacies_CamelControls = seq(0, 0.5, 0.1), ## 0 here refers to no camels, or equivalently camel interventions having zero efficacy.
		VacCampStrategies = "PROACTIVE", ImplementationDelays = 0, ImmunityDelays = 0, 
		Efficacies_Start = seq(0.00, 1, by = 0.05), ## 0 here refers to no humans, or equivalently human vaccines having zero efficacy.
		VaccineDurations = c(0, 20, 15, 10, 5, 2, 1), TimesSinceVaccination = c(0.5, 1:10))

dim(ModelRuns)
str(ModelRuns)

WriteSingleBatchFiles(ManyANDSingleFileName = ManySingleFileName)
MR_index = 1
for (MR_index in 1:dim(ModelRuns)[1]) WriteParamFile(ModelRun = ModelRuns[MR_index,])
WriteMultiBatchFile(ParamFileNames)

warnings()



## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ ## ** -- ^^ 
## ** -- ^^ REACTIVE

ManySingleFileName = "Reactive_WithCamels"	
ModelRuns = DefineModelRuns(
		VacCampStrategies = "REACTIVE", 
		Efficacies_CamelControls = seq(0.1, 0.5, 0.1),
		ImplementationDelays = seq(0,28,2), 
		ImmunityDelays = 14, 
		Efficacies_Start = seq(0.05, 1, by = 0.05), 
		ReactLevels	= c("HOSPITAL", "REGIONAL", "NATIONAL")) # c("HOSPITAL", "REGIONAL", "NATIONAL")

dim(ModelRuns)
str(ModelRuns)

WriteSingleBatchFiles(ManyANDSingleFileName = ManySingleFileName)
MR_index = 1
for (MR_index in 1:dim(ModelRuns)[1]) WriteParamFile(ModelRun = ModelRuns[MR_index,])
WriteMultiBatchFile(ParamFileNames)

warnings()








