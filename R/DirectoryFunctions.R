
PrintToDesiredPrecision = function(Vec, ndigits = 2) ### e.g. if Vec = 0.2 and ndigits = 3, will print "0.200"
{
	sprintf(paste0("%.", ndigits, "f"), Vec)
}

ChooseOutputString = function(MR, StringPrefix = "", Folder = FALSE, n_digits = 2)
{
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
			OutputString = paste0(OutputString, "_Dur_", PrintToDesiredPrecision(MR$VaccineDuration		, ndigits = n_digits))
			OutputString = paste0(OutputString, "_Lag_", PrintToDesiredPrecision(MR$TimeSinceVaccination, ndigits = n_digits))
			
		} else if (MR$VacCampStrategy == "REACTIVE")
		{
			if (MR$ReactLevel != "HOSPITAL")
			{
				if (MR$ReactLevel == "REGIONAL") OutputString = paste0(OutputString, "_reg") else
				if (MR$ReactLevel == "NATIONAL") OutputString = paste0(OutputString, "_nat")
			}
		}
		
		if (MR$Coverage				!= 1.0	)	OutputString = paste0(OutputString, "_Cov"		, PrintToDesiredPrecision(MR$Coverage, ndigits = n_digits))
		if (MR$ImplementationDelay	!= 0	)	OutputString = paste0(OutputString, "_ImpDelay"	, MR$ImplementationDelay)
		if (MR$ImmunityDelay		!= 0	)	OutputString = paste0(OutputString, "_VacDelay"	, MR$ImmunityDelay		)
		if (MR$VaccinateAllHumans			)	OutputString = paste0(OutputString, "_Blanket"							) else
		if (MR$Vaccinate_HCW				)	OutputString = paste0(OutputString, "_vHCW"								)
	}
	
	if (MR$DataFilename != "MERS_forCpp.txt")
	{
		if (MR$DataFilename == "MERS_forCpp_JanJun2013.txt"	)	OutputString = paste0(OutputString, "_JanJun2013.txt"	) else
		if (MR$DataFilename == "MERS_forCpp_JulDec2013.txt"	)	OutputString = paste0(OutputString, "_JulDec2013.txt"	) else
		if (MR$DataFilename == "MERS_forCpp_2014.txt"		)	OutputString = paste0(OutputString, "_2014.txt"			) else
																OutputString = paste0(OutputString, MR$DataFilename		)
	}
	
	if (Folder) while (substr(OutputString, 1, 1) == "_") OutputString = sub("_", "", OutputString) ### i.e. don't want underscore as first character if doing the name of a folder, so remove it. 

	return(OutputString)
}


DefineModelRuns = function( 
		
		RemoveDuplicates 			= TRUE											, 
		
		VacCampStrategies						= "REACTIVE"							,
		ReactLevels								= c("HOSPITAL")							, #c("HOSPITAL", "REGIONAL", "NATIONAL")
		Efficacies_Start						= seq(0.20, 1, by = 0.20)				,
		
		VaccineDurations						= 20000.0						, # In years. VaccineDuration == 0 => no waning. 
		TimesSinceVaccination					= 0								, # In years. 
		
		# camels
		Efficacies_CamelControls				= 0.0							,

		DataFilenames							= "MERS_forCpp.txt"				,
		
		Coverages								= 1								,
		ImplementationDelays					= 0								,
		ImmunityDelays							= 0								,
		VaccinateAllHumans_OrNot				= 0								,
		Vaccinate_HCW_OrNot						= 1								)
{
	
	### Make all possible versions of ModelRuns - note that expand.grid will not account for duplicates or for model runs that have not actually been run. 
	ModelRuns = expand.grid		(
			
			VacCampStrategy								= VacCampStrategies					,
			ReactLevel									= ReactLevels						,
			Efficacy_Start								= Efficacies_Start					,
			
			VaccineDuration								= VaccineDurations			    	, 
			TimeSinceVaccination						= TimesSinceVaccination		    	, 
			
			# camels
			Efficacy_CamelControls						= Efficacies_CamelControls			,
			
			DataFilename								= DataFilenames						,
			
			Coverage									= Coverages							,
			ImplementationDelay							= ImplementationDelays				,
			ImmunityDelay								= ImmunityDelays					,
			VaccinateAllHumans							= VaccinateAllHumans_OrNot			,
			Vaccinate_HCW								= Vaccinate_HCW_OrNot			
	)
	
	# remove scenarios where neither humans nor camels included
	ModelRuns = ModelRuns[!(ModelRuns$Efficacy_Start == 0 & ModelRuns$Efficacy_CamelControls == 0), ]
	
	# remove undefined combinations (ReactLevel irrelevant for PROACTIVE)
	ModelRuns = ModelRuns[!(ModelRuns$VacCampStrategy == "PROACTIVE" & ModelRuns$ReactLevel != "HOSPITAL"), ]
	
	# remove undefined combinations (ImplementationDelay irrelevant for PROACTIVE)
	ModelRuns = ModelRuns[!(ModelRuns$VacCampStrategy == "PROACTIVE" & ModelRuns$ImplementationDelay != 0), ]
	
	#### Create the OutputStrings for all possible model runs...
	OutputFolderNames 	= rep(NA, dim(ModelRuns)[1])
	for (run in 1:dim(ModelRuns)[1])	OutputFolderNames[run] = ChooseOutputString(ModelRuns[run,], Folder = TRUE)
	#### .... bind them to data.frame
	ModelRuns = cbind(ModelRuns, OutputFolderNames)
	
	#### Remove duplicates
	if (RemoveDuplicates)			ModelRuns = ModelRuns[!duplicated(ModelRuns$OutputFolderNames), ]
	
	return(ModelRuns)
}





