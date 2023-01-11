MAX_ONSET_DAY = 549

##### ===== ..... ##### ===== ..... ##### ===== ..... ##### ===== ..... ##### ===== ..... ##### ===== ..... ##### ===== ..... ##### ===== ..... ##### ===== ..... ##### ===== .....  
##### after everything put in a data frame

PrintToDesiredPrecision = function(Vec, ndigits = 2) ### e.g. if Vec = 0.2 and ndigits = 3, will print "0.200"
{
	sprintf(paste0("%.", ndigits, "f"), Vec)
}
GetMinDate = function(MR_index)
{
	return(min(ModelRuns[MR_index,]$Pruning_Range[[1]]))
}
GetMaxDate = function(MR_index)
{
	return(max(ModelRuns[MR_index,]$Pruning_Range[[1]]))
}

ChooseOutputString = function(MR_index, StringPrefix = "", Folder = FALSE, n_digits = 2)
{
	MR = ModelRuns[MR_index, ]
	
	if (!(MR$VacCampStrategy %in% c("REACTIVE", "PROACTIVE")))
		stop("ChooseOutputString error: MR$VacCampStrategy not reactive")
	
	OutputString = StringPrefix
	
	if (MR$Efficacy_Start == 0) OutputString = paste0(OutputString, "_NoHumans") else
								OutputString = paste0(OutputString, "_Eff_", PrintToDesiredPrecision(MR$Efficacy_Start, ndigits = n_digits))
	
	if (MR$Efficacy_CamelControls > 0)
		OutputString = paste0(OutputString, "_CamelsControlEff_", PrintToDesiredPrecision(MR$Efficacy_CamelControls, ndigits = n_digits))
	
	if (MR$Efficacy_Start != 0) ##  i.e. if vaccinating humans
	{
		if (MR$VacCampStrategy == "PROACTIVE")
		{
			OutputString = paste0(OutputString, "_ProAct")
			if (!MR$ExpWaning)	OutputString = paste0(OutputString, "_Hill_", PrintToDesiredPrecision(MR$HillPower, 0))
			
			OutputString = paste0(OutputString, "_Dur_", PrintToDesiredPrecision(MR$VaccineDuration		, ndigits = n_digits))
			OutputString = paste0(OutputString, "_Lag_", PrintToDesiredPrecision(MR$TimeSinceVaccination, ndigits = n_digits))
			
		} else if (MR$VacCampStrategy == "REACTIVE")
		{
			if (MR$WaningInReactive & MR$VaccineDuration != 0) 
			{
				OutputString = paste0(OutputString, "_RW_Dur_", PrintToDesiredPrecision(MR$VaccineDuration, 2))
				if (!MR$ExpWaning)	OutputString = paste0(OutputString, "_Hill_", PrintToDesiredPrecision(MR$HillPower, 0))
			}
			if (MR$ReactLevel != "HOSPITAL")
			{
				if (MR$ReactLevel == "REGIONAL") OutputString = paste0(OutputString, "_reg") else
				if (MR$ReactLevel == "NATIONAL") OutputString = paste0(OutputString, "_nat")
			}
			if (MR$DoTriggers) 
				OutputString = paste0(OutputString, 
						"_Tr", PrintToDesiredPrecision(MR$TrigThreshold, 0), 
						"_",   PrintToDesiredPrecision(MR$TrigTimeframe, 0))
		}
		
		if (MR$Coverage				!= 1.0	)	OutputString = paste0(OutputString, "_Cov"		, PrintToDesiredPrecision(MR$Coverage, ndigits = n_digits))
		if (MR$ImplementationDelay	!= 0	)	OutputString = paste0(OutputString, "_ImpDelay"	, MR$ImplementationDelay)
		if (MR$ImmunityDelay		!= 0	)	OutputString = paste0(OutputString, "_VacDelay"	, MR$ImmunityDelay		)
		if (MR$VaccinateAllHumans			)	OutputString = paste0(OutputString, "_Blanket"							) else
		if (MR$Vaccinate_HCW				)	OutputString = paste0(OutputString, "_vHCW"								)
	}
	
	MinDay_Pruning = min(MR$Pruning_Range[[1]])
	MaxDay_Pruning = max(MR$Pruning_Range[[1]])
	
	if (MinDay_Pruning > 0 | MaxDay_Pruning < (MAX_ONSET_DAY - 1))
	{
		OutputString = paste0(OutputString, "_s", MinDay_Pruning) # start date
		OutputString = paste0(OutputString, "_f", MaxDay_Pruning) # end date
	}
	
	if (Folder) while (substr(OutputString, 1, 1) == "_") OutputString = sub("_", "", OutputString) ### i.e. don't want underscore as first character if doing the name of a folder, so remove it. 

	return(OutputString)
}


DefineModelRuns = function( 
		
		RemoveDuplicates 			= TRUE											, 
		IncludeCompletedRunsOnly 	= TRUE						, 
		
		VacCampStrategies			= "REACTIVE"							,
		WaningInReactiveAndNot		= 0										,
		ReactLevels					= c("HOSPITAL")							, #c("HOSPITAL", "REGIONAL", "NATIONAL")
		Efficacies_Start			= seq(0.20, 1, by = 0.20)				,
		
		VaccineDurations			= 20000.0						, # In years. VaccineDuration == 0 => no waning. 
		TimesSinceVaccination		= 0								, # In years.
		
		ExpWaningAndNot				= 1								,
		HillPowers					= 4								,
		
		Pruning_Ranges				= list( c(0, MAX_ONSET_DAY - 1)),
		
		# camels
		Efficacies_CamelControls	= 0.0							,
		
		## Triggers
		DoTriggersAndNot			= 0	, 
		TrigThresholds				= 1	,
		TrigTimeframes				= 1	,
		
		Coverages					= 1								,
		ImplementationDelays		= 0								,
		ImmunityDelays				= 0								,
		VaccinateAllHumans_OrNot	= 0								,
		Vaccinate_HCW_OrNot			= 1								)
{
	### Make all possible versions of ModelRuns - note that expand.grid will not account for duplicates or for model runs that have not actually been run. 
	ModelRuns = expand.grid		(
			
			VacCampStrategy			= VacCampStrategies		,
			ReactLevel				= ReactLevels			,
			WaningInReactive		= WaningInReactiveAndNot,
			Efficacy_Start			= Efficacies_Start		,
			
			VaccineDuration			= VaccineDurations			, 
			TimeSinceVaccination	= TimesSinceVaccination		, 
			
			ExpWaning				= ExpWaningAndNot		,
			HillPower				= HillPowers				,
			
			Pruning_Range			= Pruning_Ranges		,
			
			# camels
			Efficacy_CamelControls	= Efficacies_CamelControls,
			
			## Triggers
			DoTriggers				= DoTriggersAndNot		, 
			TrigThreshold			= TrigThresholds		,
			TrigTimeframe 			= TrigTimeframes		,
			
			Coverage				= Coverages				,
			ImplementationDelay		= ImplementationDelays	,
			ImmunityDelay			= ImmunityDelays		,
			VaccinateAllHumans		= VaccinateAllHumans_OrNot,
			Vaccinate_HCW			= Vaccinate_HCW_OrNot			
	)
	cat("expanded grid, ")
	
	# remove scenarios where neither humans nor camels included
	ModelRuns = ModelRuns[!(ModelRuns$Efficacy_Start == 0 & ModelRuns$Efficacy_CamelControls == 0), ]
	
	# remove undefined combinations (ReactLevel irrelevant for PROACTIVE)
	ModelRuns = ModelRuns[!(ModelRuns$VacCampStrategy == "PROACTIVE" & ModelRuns$ReactLevel != "HOSPITAL"), ]
	
	# remove undefined combinations (WaningInReactive irrelevant for PROACTIVE)
	ModelRuns = ModelRuns[!(ModelRuns$VacCampStrategy == "PROACTIVE" & ModelRuns$WaningInReactive != 0), ]
	
	# remove undefined combinations (ImplementationDelay irrelevant for PROACTIVE)
	ModelRuns = ModelRuns[!(ModelRuns$VacCampStrategy == "PROACTIVE" & ModelRuns$ImplementationDelay != 0), ]
	
	# remove undefined combinations (ImmunityDelay == 14 for REACTIVE)
	ModelRuns = ModelRuns[!(ModelRuns$VacCampStrategy == "REACTIVE" & ModelRuns$ImmunityDelay == 0), ]
	
	cat("undefined scenarios removed, ")
	if (dim(ModelRuns)[1] == 0) stop("ModelRuns empty")
	
	ModelRuns <<- ModelRuns
	
	#### Create  OutputStrings for all possible model runs...
	cat("Creating output strings, ")
	OutputFolderNames 	= mapply(ChooseOutputString	, MR_index = 1:dim(ModelRuns)[1], MoreArgs = list(Folder = TRUE))
	cat("choosing min days, ")
	MinDay_Pruning		= mapply(GetMinDate			, MR_index = 1:dim(ModelRuns)[1])
	cat("choosing max days, ")
	MaxDay_Pruning		= mapply(GetMaxDate			, MR_index = 1:dim(ModelRuns)[1])
	cat("loop finished, ")
	
	#### .... bind them to data.frame
	ModelRuns = cbind(ModelRuns, MinDay_Pruning, MaxDay_Pruning, OutputFolderNames)
	
	#### Remove duplicates
	if (RemoveDuplicates)			ModelRuns = ModelRuns[!duplicated(ModelRuns$OutputFolderNames), ]
	#### Remove incomplete runs, if you want (e.g. useful to not compare runs you haven't done)
	
	# Check if either full chains or summary chains available
	RunsCompleted = file.exists(file.path(CppOutputDirectory, paste0("CF_Chains_"		, ModelRuns$OutputFolderNames, ".txt"))) 	| 
					file.exists(file.path(CppOutputDirectory, paste0("CF_ChainsSumm_"	, ModelRuns$OutputFolderNames, ".txt")))
	if (IncludeCompletedRunsOnly) 	ModelRuns = ModelRuns[RunsCompleted, ] ### _underscore after ParameterChainOutput depends on Folder argument in ChooseOutputString function called above.
	
	cat("done.")
	
	return(ModelRuns)
}





