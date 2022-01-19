
#include "InputOutput.h"
#include "Globals.h"
#include "Likelihood.h"
#include "MCMC.h"
#include "DJL_Structs.h"
#include "Macros.h"

VaccCampaignStrategy	Convert_VaccCampaignStrategy_FromString	(const std::string& OptionString)
{
		 if (OptionString == "REACTIVE"		) return VaccCampaignStrategy::REACTIVE;
	else if (OptionString == "PROACTIVE"	) return VaccCampaignStrategy::PROACTIVE;
	else 
	{
		std::cout << endl << "Convert_VaccCampaignStrategy_FromString ERROR: String not recognized" << endl;
		exit(1);
	}
}
ReactiveLevel			Convert_ReactiveLevel_FromString		(const std::string& OptionString)
{
		 if (OptionString == "HOSPITAL") return ReactiveLevel::HOSPITAL;
	else if (OptionString == "REGIONAL") return ReactiveLevel::REGIONAL;
	else if (OptionString == "NATIONAL") return ReactiveLevel::NATIONAL;
	else 
	{
		std::cout << endl << "Convert_ReactiveLevel_FromString ERROR: String not recognized" << endl;
		exit(1);
	}
}
std::string Convert_VaccCampaignStrategy_FromEnumClass	(VaccCampaignStrategy VacCampStrategy)
{
		 if (VacCampStrategy == VaccCampaignStrategy::REACTIVE	) return "REACTIVE"		;
	else if (VacCampStrategy == VaccCampaignStrategy::PROACTIVE	) return "PROACTIVE"	;	
	else 
	{
		std::cout << endl << "Convert_VaccCampaignStrategy_FromString ERROR: Class not recognized" << endl;
		exit(1);
	}
}
std::string Convert_ReactiveLevel_FromEnumClass			(ReactiveLevel ReactLevel)
{
		 if (ReactLevel == ReactiveLevel::HOSPITAL) return "HOSPITAL";
	else if (ReactLevel == ReactiveLevel::REGIONAL) return "REGIONAL";
	else if (ReactLevel == ReactiveLevel::NATIONAL) return "NATIONAL";
	else 
	{
		std::cout << endl << "Convert_VaccCampaignStrategy_FromEnumClass ERROR: Class not recognized" << endl;
		exit(1);
	}
}

bool FileExists				(string FileName)
{
	std::ifstream infile(FileName);
	bool FileOkay = infile.good();
	infile.close(); 
	return FileOkay;
}
void ReadInParams			(ModelRun &MR, string pParamFileName)
{
	ifstream ParamsEtc; 
	std::cout << "Reading Param file: " << pParamFileName ; 
	fflush(stderr);	fflush(stdout);

	string param_name, param_value_string;
	if (FileExists(pParamFileName))
	{
		ParamsEtc.open(pParamFileName);
		while (!ParamsEtc.eof())
		{
			std::getline(ParamsEtc, param_name			, '\t');
			std::getline(ParamsEtc, param_value_string	, '\n');

				 if (param_name == "NumIterations"					)	MR.NumIterations				= std::stoi(param_value_string);	// int / bool
			else if (param_name == "BurnIn"							)	MR.BurnIn						= std::stoi(param_value_string);	// int / bool
			else if (param_name == "StoreEvery"						)	MR.StoreEvery					= std::stoi(param_value_string);	// int / bool
			else if (param_name == "DJL_InputData"					)	MR.DJL_InputData				= std::stoi(param_value_string);	// int / bool
			else if (param_name == "OutputTreesEvery"				)	MR.OutputTreesEvery				= std::stoi(param_value_string);	// int / bool
			else if (param_name == "Efficacy_Start"					)	MR.Efficacy_Start				= std::stod(param_value_string);	// double
			else if (param_name == "VaccineDuration"				)	MR.VaccineDuration				= std::stod(param_value_string);	// double
			else if (param_name == "TimeSinceVaccination"			)	MR.TimeSinceVaccination			= std::stod(param_value_string);	// double
			else if (param_name == "Coverage"						)	MR.Coverage						= std::stod(param_value_string);	// double
			else if (param_name == "ImplementationDelay"			)	MR.ImplementationDelay			= std::stoi(param_value_string);	// int / bool
			else if (param_name == "ImmunityDelay"					)	MR.ImmunityDelay				= std::stoi(param_value_string);	// int / bool
			else if (param_name == "VaccinateAllHumans"				)	MR.VaccinateAllHumans			= std::stoi(param_value_string);	// int / bool
			else if (param_name == "Vaccinate_HCW"					)	MR.Vaccinate_HCW				= std::stoi(param_value_string);	// int / bool
			else if (param_name == "Efficacy_CamelControls"			)	MR.Efficacy_CamelControls		= std::stod(param_value_string);	// double
			else if (param_name == "VacCampStrategy"				)	MR.VacCampStrategy				= Convert_VaccCampaignStrategy_FromString	(param_value_string);
			else if (param_name == "ReactLevel"						)	MR.ReactLevel					= Convert_ReactiveLevel_FromString			(param_value_string);
			else	std::cout << " ReadInParams ERROR: param_name " << param_name << " not recognized. param_value_string = " << param_value_string << endl;
		}
		ParamsEtc.close();
	}
	std::cout << " DONE" << std::endl;
	fflush(stderr);	fflush(stdout);
}
//// Colname functions
void PutColNamesOn_Chains		(ofstream &Chains)
{
	Chains << "LogLikelihood"					<< "\t";		// Log likelihood
	Chains << "Mean_GenTime"					<< "\t";		// Mean generation time (taken from Lancet ID)
	Chains << "SD_GenTime"						<< "\t";		// Generation time std (taken from Lancet ID)
	Chains << "R_Hosp"							<< "\t";		// R hosp
	Chains << "R_SameRegion"					<< "\t";		// R region
	Chains << "R_OtherRegion"					<< "\t";		// R other region
	Chains << "k_Hosp"							<< "\t";		// k hosp  //_cvParam = 0: param5 =  k ; _cvParam = 1: param5 = cv = 1/sqrt(k)
	Chains << "k_Region"						<< "\t";		// k region
	Chains << "Mean_NumIntro_atZero"			<< "\t";		// Mean number of introductions at time 0
	Chains << "Mean_NumIntro_ExpRate"			<< "\t";		// Mean number of introductions - exponential rate
	Chains << "Overdispersion_Intro"			<< "\t";		// overdispersion intro
	Chains << "Reduction_ClusterTrans"			<< "\t";		// reduction in cluster transmission over time
	Chains << "PowerReduction_ClusterTrans"		<< "\t";		// power reduction in cluster transmission over time
	Chains << "PeakSeasonalityTime"				<< "\t";		// Time of year when seasonality peak (0 = 1st Jan), (1 = 31st Dec)
	Chains << "DK_Param_1"						<< "\t";
	Chains << "DK_Param_2"						<< "\t";
	Chains << "ProbSource_Hosp"					<< "\t";		
	Chains << "ProbSource_SameRegion"			<< "\t";		
	Chains << "ProbSource_OtherRegion"			<< "\t";		
	Chains << "ProbSource_Reservoir"			<< "\t";		
	Chains << "CumNumIntroducers"				<< "\t";		// Total number of introducers summed over all days. 
	Chains << "_simulNbIntro"					<< "\t";
	Chains << "_simulNbCase"					<< "\t";  
	Chains << "_simulNbCluster"					<< "\t"; 
	Chains << "_simulMeanClusterSize"			<< "\t"; 
	Chains << "_simulMaxClusterSize"			<< "\t";
	Chains << "_simulPClusterSize1"				<< "\t";
	Chains << "_simulPClusterSize10"			<< "\t";
	Chains << "_simulMaxNbCaseOverPeriod"		<< "\t";
	Chains << "_simulMaxNbClusterOverPeriod"	<< std::endl;
}
void PutColNamesOn_IndividualR	(ofstream &IndividualR)
{
	IndividualR << "CaseNumber"				<< "\t"; 
	IndividualR << "region"					<< "\t"; 
	IndividualR << "hospital"				<< "\t"; 
	IndividualR << "clusterOfCase"			<< "\t"; 
	IndividualR << "onset"					<< "\t"; 
	IndividualR << "rankOfCaseInCluster"	<< "\t";
	IndividualR << "IndObsR_Hosp"			<< "\t";
	IndividualR << "IndObsR_SameRegion"		<< "\t";
	IndividualR << "IndObsR_OtherRegion"	<< "\t";
	IndividualR << endl;
}
void PutColNamesOn_Trees		(ofstream &Trees)
{
	Trees << "iteration"		<< "\t";
	Trees << "CaseNumber"		<< "\t";
	Trees << "hospital"			<< "\t"; 
	Trees << "region"			<< "\t"; 
	Trees << "onset"			<< "\t"; 
	Trees << "H2Htype"			<< "\t"; 
	Trees << "Infector"			<< "\t";
	Trees << "hospital_Infector"<< "\t"; 
	Trees << "Onset_Infector"	<< "\t"; 
	Trees << endl;
}
void PutColNamesOn_Trees_CF		(ofstream &outputFileTree_CF)
{
	outputFileTree_CF << "iteration"			<< "\t";
	outputFileTree_CF << "CFnum"				<< "\t";
	outputFileTree_CF << "CaseNumber"			<< "\t";
	outputFileTree_CF << "hospital"				<< "\t"; 
	outputFileTree_CF << "onset"				<< "\t"; 
	outputFileTree_CF << "H2Htype"				<< "\t"; 
	outputFileTree_CF << "Infector"				<< "\t";
	outputFileTree_CF << "hospital_Infector"	<< "\t"; 
	outputFileTree_CF << "Onset_Infector"		<< "\t"; 
	outputFileTree_CF << endl;
}
void PutColNamesOn_SimulCluster	(ofstream &outputFileSimulCluster)
{
	////// set column names for outputFileSimulCluster 
	outputFileSimulCluster << "iteration" << "\t" << "_simulClusterStart" << "\t" << "_regionOfCluster" << "\t" << "_hospitalOfCluster" << "\t" << "_casesInCluster" << "\t" << "_simulClusterR" << endl;
}
void PutColNamesOn_AllOutput	(AllOutput& OUTPUT)
{
	if (OUTPUT.Write_Chains			)	PutColNamesOn_Chains		(OUTPUT.Chains			);
	if (OUTPUT.Write_IndividualR	)	PutColNamesOn_IndividualR	(OUTPUT.IndividualR		);
	if (OUTPUT.Write_Trees			)	PutColNamesOn_Trees			(OUTPUT.Trees			);
	if (OUTPUT.Write_Trees_CF		)	PutColNamesOn_Trees_CF		(OUTPUT.Trees_CF		);
	if (OUTPUT.Write_SimulCluster	)	PutColNamesOn_SimulCluster	(OUTPUT.SimulCluster	);
}	
void InitializeOutputValues		(AllOutput& OUTPUT)
{
	//// must be called after loadAndInitializeData function has been called. 
	if (OUTPUT.Write_ClusterR)
	{
		for (int iCluster = 0; iCluster < _numberOfClusters; iCluster++) OUTPUT.ClusterR << _hospitalOfCluster[iCluster]					<< "\t";	OUTPUT.ClusterR << endl;
		for (int iCluster = 0; iCluster < _numberOfClusters; iCluster++) OUTPUT.ClusterR << _regionOfHospital[_hospitalOfCluster[iCluster]]	<< "\t";	OUTPUT.ClusterR << endl;
		for (int iCluster = 0; iCluster < _numberOfClusters; iCluster++) OUTPUT.ClusterR << _startOfCluster[iCluster]						<< "\t";	OUTPUT.ClusterR << endl;
	}
	if (OUTPUT.Write_ClusterReductionR)
	{
		for (int iCluster = 0; iCluster < _numberOfClusters; iCluster++) OUTPUT.ClusterReductionR	<< _hospitalOfCluster[iCluster]						<< "\t";	OUTPUT.ClusterReductionR << endl;
		for (int iCluster = 0; iCluster < _numberOfClusters; iCluster++) OUTPUT.ClusterReductionR	<< _regionOfHospital[_hospitalOfCluster[iCluster]]	<< "\t";	OUTPUT.ClusterReductionR << endl;
		for (int iCluster = 0; iCluster < _numberOfClusters; iCluster++) OUTPUT.ClusterReductionR	<< _startOfCluster[iCluster]						<< "\t";	OUTPUT.ClusterReductionR << endl;
	}

	if (OUTPUT.Write_SimulCluster)
		for (int iCluster = 0; iCluster < _numberOfClusters; iCluster++)
		{
			OUTPUT.SimulCluster << 0 << "\t" << _startOfCluster[iCluster] << "\t" << _regionOfHospital[_hospitalOfCluster[iCluster]] << "\t";
			OUTPUT.SimulCluster << _hospitalOfCluster[iCluster] << "\t" << _casesInCluster[iCluster].size() << "\t" << 0 << endl;
		}
}

void loadAndInitializeData	(string fileName, ModelRun& MR, bool printData = false)
{
	std::cout << "loadAndInitializeData: " << std::endl;
	fflush(stderr);	fflush(stdout);

	// NOTE: _maxDuration, _beta, _parameter[0], and _parameter[1] must be set BEFORE calling this function
	int lineCount = 0, dummy;
	char line[2000];

	// Check that file is there
	ifstream iFile(fileName.c_str());
	if (!iFile.is_open())
	{
		cout << "ERROR: file " << fileName << " does not exist!" << endl;
		exit(1);
	}
	// Determine number of lines
	while (!iFile.eof())
	{
		iFile.getline(line, 2000);
		++lineCount;
	}
	cout << "Number of lines read: " << lineCount << endl;
	fflush(stderr);
	fflush(stdout);

	// Rewind to beginning of the file
	iFile.clear(); // Clear fail and eof bits
	iFile.seekg(0, ios_base::beg); // Back to the start

	// Initialize data
	_numberOfCases			= lineCount + 1; // Community is case 0
	if (MR.DJL_InputData)	_numberOfCases = lineCount - 1;

	_nameCase				= vector<int>(_numberOfCases, 0);
	_onset					= vector<int>(_numberOfCases, 0);
	_onset_week				= vector<int>(_numberOfCases, 0);
	_hospital				= vector<int>(_numberOfCases, 0);
	_region					= vector<int>(_numberOfCases, 0);
	_infector				= vector<int>(_numberOfCases, 0);
	_individualObsR			= vector<vector<int>>	(_numberOfCases, vector<int>	(3, 0));
	_meanIndividualObsR		= vector<vector<double>>(_numberOfCases, vector<double>	(3, 0));
	_individualExpR			= vector<vector<double>>(_numberOfCases, vector<double> (3, 1));
	_numberOfIntroducers	= map<int, int>();
	_generation				= vector<int>(_numberOfCases, -1);
	_secondaryCases			= vector<set<int>>(_numberOfCases, set<int>());

	//// DJL containters
	_vaccinated				= vector<bool>(_numberOfCases, 0);
	_protected				= vector<bool>(_numberOfCases, 0);
	_DeleteCase_CF			= vector<bool>(_numberOfCases, 0);
	_HealthCareWorker		= vector<int> (_numberOfCases, 0);
	_Sex					= vector<int> (_numberOfCases, 0); 
	_Dead					= vector<int> (_numberOfCases, 0); 
	_Symptomatic			= vector<int> (_numberOfCases, 0); 
	_Age					= vector<int> (_numberOfCases, 0); 

	// Read again to load data
	if (MR.DJL_InputData)
	{
		int NumDataVariables = 10; 
		std::vector<string> DataNames(NumDataVariables, "");
		for (int VarName = 0; VarName < NumDataVariables; VarName++) iFile >> DataNames[VarName];
		for (int VarName = 0; VarName < NumDataVariables; VarName++) std::cout << " " << DataNames[VarName]; 
		std::cout << std::endl; 
	}
	if (MR.DJL_InputData)
	{
		for (int iCase = 1; iCase < _numberOfCases; ++iCase)
		{
			iFile >> _nameCase	[iCase];
			iFile >> _onset		[iCase];
			iFile >> _hospital	[iCase];			// hospital
			iFile >> _region	[iCase];			// region
			iFile >> dummy;							// city			
			iFile >> _HealthCareWorker[iCase];		// HCW
			iFile >> _Sex[iCase];					// Sex
			iFile >> _Dead[iCase];					// Dead
			iFile >> _Symptomatic[iCase];			// Symptomatic
			iFile >> _Age[iCase];					// Age in years
		}
	}
	else
	{
		for (int iCase = 1; iCase < _numberOfCases; ++iCase)
		{
			iFile >> _nameCase[iCase];
			iFile >> _onset[iCase];
			iFile >> _hospital[iCase];		// cluster
			iFile >> dummy;					// city			//// hospital.city in line list
			iFile >> _region[iCase];		// region
		}
	}
	for (int iCase = 1; iCase < _numberOfCases; ++iCase) _onset_week[iCase] = _onset[iCase] / 7; 

	iFile.close();

	// Find first cases in hospital, region and country
	_FirstOnsetInHosp = vector<int>(_numberOfHospitals, 100000); //// i.e. initialize vector with arbitrarily large number
	_FirstOnsetInRegion = vector<int>(_numberOfRegions, 100000); //// i.e. initialize vector with arbitrarily large number
	_FirstOnsetInCountry = 100000; 
	_numberOfHCWs = 0;
	for (int iCase = 1; iCase < _numberOfCases; ++iCase)
	{
		if (_onset[iCase] < _FirstOnsetInHosp[_hospital[iCase]]	)	_FirstOnsetInHosp[_hospital[iCase]] = _onset[iCase]; //// Find first case in each hospital. 
		if (_onset[iCase] < _FirstOnsetInRegion[_region[iCase]]	)	_FirstOnsetInRegion[_region[iCase]] = _onset[iCase]; //// Find first case in each region. 
		if (_onset[iCase] < _FirstOnsetInCountry				)	_FirstOnsetInCountry				= _onset[iCase]; //// Find first case in each region. 
		if (_HealthCareWorker[iCase] == 1) _numberOfHCWs++; 
	}

	_minDay = 0;	//*min_element(_onset.begin(), _onset.end());
	_maxDay = *max_element(_onset.begin(), _onset.end());

	// Initialize _possibleInfectors (0 = Reservoir is always a possible infector). Possible infectors are set here then never altered throughout program. 
	_possibleInfectors = vector<set<int>>(_numberOfCases, set<int>({0}));
	bool check1, check2;
	for (int iCase = 1; iCase < _numberOfCases; ++iCase)
		for (int currCase = 1; currCase < _numberOfCases; ++currCase)
		{
			if (currCase == iCase) continue;
			check1 = _onset[iCase] > _onset[currCase];								//// was infectee (iCase) infected AFTER currCase (potential infector)? 
			check2 = _onset[iCase] - _onset[currCase] < _maxDuration;				//// was time difference between onset for infectee (iCase) and onset for currCase (potential infector) less than maximum possible difference? 
			if (check1 && check2) _possibleInfectors[iCase].insert(currCase);		//// if yes to both questions, add to list of possible infectors for this infectee. 
		}

	// Initialize _numberOfIntroducers
	for (int day = _minDay; day <= _maxDay; ++day)	_numberOfIntroducers[day] = 0;
	_numberOfIntroducersWith0Sec = 0;

	// set _regionOfHospital
	_regionOfHospital = vector<int>(_numberOfHospitals, -1);
	for (int iCase = 1; iCase<_numberOfCases; iCase++)
		_regionOfHospital[_hospital[iCase]] = _region[iCase];

	_numberOfHospitalsPerRegion = vector<int>(_numberOfRegions, 0);
	_hospitalInRegion			= vector<vector<int>>(_numberOfRegions, vector<int>());
	for (int hospital = 0; hospital < _numberOfHospitals; hospital++)
	{
		_numberOfHospitalsPerRegion[_regionOfHospital[hospital]]++;
		_hospitalInRegion[_regionOfHospital[hospital]].push_back(hospital);
	}

	_minDelayBetweenClusters	= _maxDuration;
	_clusterOfCase				= vector<int>(_numberOfCases, -1);

	int currentCluster = -1;
	_numberOfClusters = 0;
	int dayLastCase, delayToLastCase;
	for (int hospital = 0; hospital < _numberOfHospitals; hospital++)
	{
		dayLastCase = -1000;
		for (int iCase = 1; iCase < _numberOfCases; ++iCase)
			if (_hospital[iCase] == hospital)
			{
				delayToLastCase = _onset[iCase] - dayLastCase;
				if (delayToLastCase > _minDelayBetweenClusters)
				{
					currentCluster++;
					_numberOfClusters++;
					_startOfCluster.push_back(_onset[iCase]);
					_hospitalOfCluster.push_back(hospital);
				}
				_clusterOfCase[iCase] = currentCluster;
				dayLastCase = _onset[iCase];
			}
	}

	//// set cases in Region and cluster
	_casesInRegion	= vector<set<int>>(_numberOfRegions	, set<int>({}));
	_casesInCluster = vector<set<int>>(_numberOfClusters, set<int>({}));
	for (int iCase = 1; iCase < _numberOfCases; ++iCase) _casesInRegion	[_region		[iCase]].insert(iCase);
	for (int iCase = 1; iCase < _numberOfCases; ++iCase) _casesInCluster[_clusterOfCase	[iCase]].insert(iCase);

	_rankOfCaseInCluster = vector<int>(_numberOfCases, -1);	// more specifically: nb of cases in cluster prior to onset of this case
	for (int iCase = 1; iCase < _numberOfCases; ++iCase)
	{
		int currentRank = 0;
		for (set<int>::iterator it = _casesInCluster[_clusterOfCase[iCase]].begin(); it != _casesInCluster[_clusterOfCase[iCase]].end(); ++it)
			if (_onset[*it] >= _onset[iCase]) break;		else currentRank++;
		_rankOfCaseInCluster[iCase] = currentRank;
	}

	_clusterExpR		= vector<vector<double>>(_numberOfClusters	, vector<double>(3, 1)); //// vector of vector of doubles, however second index only ever takes value 0 = SameHosptial. 
	_clusterReductionR	= vector<double>		(_numberOfClusters	, 1);
	_regionExpR			= vector<double>		(_numberOfRegions	, 1);


	// Initialize GT density and weighted proba
	_densityGT					= vector<double>		(_maxDuration	, 0.);
	_probaOfTimingTrans			= vector<vector<double>>(_numberOfCases	, vector<double>(_numberOfCases, 0.)); //// i.e. an _numberOfCases by _numberOfCases matrix initialized to zero. 
	_weightedProbaOfTimingTrans = vector<vector<double>>(_numberOfCases	, vector<double>(_numberOfCases, 0.)); //// i.e. an _numberOfCases by _numberOfCases matrix initialized to zero. 
	updateDensityGT();
	updateProbaOfTimingTrans();
	// Assign missing infectors
	initiateInfector();
	if (_variabilityRAtClusterLevel > 0) initiateVariabilityExpR();

	// Print data
	if (printData)
	{
		cout << "-----------------------------------------------" << endl;
		for (int iCase = 0; iCase < _numberOfCases; ++iCase)
		{
			cout << _nameCase[iCase] << " onset " << _onset[iCase];
			cout << ", probably infected by " << _nameCase[_infector[iCase]] << ",";
			cout << " generation: " << _generation[iCase];
			cout << " infector of: ";
			for (set<int>::iterator it = _secondaryCases[iCase].begin(); it != _secondaryCases[iCase].end(); ++it) cout << _nameCase[*it] << " ";
			cout << endl;
		}
	}

	_simulCaseOnDay				= vector<set<int>>	(_maxDay, set<int>());
	_simulLastOnsetInHospital	= vector<int>		(_numberOfHospitals, -1000);
	_simulRegionR				= vector<double>	(_numberOfRegions);

	std::cout << "loadAndInitializeData DONE " << std::endl;
	fflush(stderr);	fflush(stdout);
}

void WriteModelMetaData(AllOutput &OUTPUT, ModelRun &MR, FileStrings_Struct FileStrings, double delta, double kRinit, double pDetectIndex, double kIntro)  //// unseen "arguments" are various global variables, e.g. _spatialLevel
{
	//// ModelRunVariables
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.NumIterations		) << "\t" << MR.NumIterations		<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.BurnIn				) << "\t" << MR.BurnIn				<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.NumCFQuantities		) << "\t" << MR.NumCFQuantities		<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.OutputTreesEvery	) << "\t" << MR.OutputTreesEvery	<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.NumCFsPerTree		) << "\t" << MR.NumCFsPerTree		<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.StoreEvery			) << "\t" << MR.StoreEvery			<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.RunOnCluster		) << "\t" << MR.RunOnCluster		<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.DJL_InputData		) << "\t" << MR.DJL_InputData		<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.max_threads			) << "\t" << MR.max_threads			<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.seed1				) << "\t" << MR.seed1				<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.seed2				) << "\t" << MR.seed2				<< std::endl;

	//// DJL ModelRun
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.Efficacy_Start				)	<< "\t" << MR.Efficacy_Start			<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.Efficacy_Current			)	<< "\t" << MR.Efficacy_Current			<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.Coverage					)	<< "\t" << MR.Coverage					<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.ImplementationDelay			)	<< "\t" << MR.ImplementationDelay		<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.VaccineDelay				)	<< "\t" << MR.ImmunityDelay				<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.VaccinateAllHumans			)	<< "\t" << MR.VaccinateAllHumans		<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.Vaccinate_HCW				)	<< "\t" << MR.Vaccinate_HCW				<< std::endl;

	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.VacCampStrategy	)	<< "\t" << Convert_VaccCampaignStrategy_FromEnumClass	(MR.VacCampStrategy	) << std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.ReactLevel		)	<< "\t" << Convert_ReactiveLevel_FromEnumClass			(MR.ReactLevel		) << std::endl;

	// camels
	OUTPUT.MetaData << GET_VARIABLE_NAME(MR.Efficacy_CamelControls		)	<< "\t" << MR.Efficacy_CamelControls		<< std::endl;
	
	///// FileStrings
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.Simon_scenarioName		)	<< "\t" << FileStrings.Simon_scenarioName	<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.DJL_scenarioName		)	<< "\t" << FileStrings.DJL_scenarioName		<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.inputFileName			)	<< "\t" << FileStrings.inputFileName		<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.inputFile				)	<< "\t" << FileStrings.inputFile			<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.Chains					)	<< "\t" << FileStrings.Chains				<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.IndividualR			)	<< "\t" << FileStrings.IndividualR			<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.Infector				)	<< "\t" << FileStrings.Infector				<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.ClusterR				)	<< "\t" << FileStrings.ClusterR				<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.ClusterReductionR		)	<< "\t" << FileStrings.ClusterReductionR	<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.SimulCluster			)	<< "\t" << FileStrings.SimulCluster			<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.Tree					)	<< "\t" << FileStrings.Tree					<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.Tree_CF				)	<< "\t" << FileStrings.Tree_CF				<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.CF_Chains				)	<< "\t" << FileStrings.CF_Chains			<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.MetaData				)	<< "\t" << FileStrings.MetaData				<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(FileStrings.InitParamValues		)	<< "\t" << FileStrings.InitParamValues		<< std::endl;
																						   		
	///// Global Variables
	OUTPUT.MetaData << GET_VARIABLE_NAME(_dataSevereCases					)	<< "\t" << _dataSevereCases							<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_withReductionR					)	<< "\t" << _withReductionR							<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_withSeasonality					)	<< "\t" << _withSeasonality							<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_introOverdispersed				)	<< "\t" << _introOverdispersed						<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_estimOverdispersionIntro			)	<< "\t" << _estimOverdispersionIntro				<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_variabilityRAtClusterLevel 		)	<< "\t" << _variabilityRAtClusterLevel 				<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_spatialLevel						)	<< "\t" << _spatialLevel							<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_heterogeneityR					)	<< "\t" << _heterogeneityR							<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_heterogeneityReductionR			)	<< "\t" << _heterogeneityReductionR					<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_maxDuration						)	<< "\t" << _maxDuration								<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_trackGeneration					)	<< "\t" << _trackGeneration							<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_kRfix								)	<< "\t" << _kRfix									<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_heterogeneityRegionR				)	<< "\t" << _heterogeneityRegionR					<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_cvParam							)	<< "\t" << _cvParam									<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_numberOfCases						)	<< "\t" << _numberOfCases					<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_numberOfHospitals					)	<< "\t" << _numberOfHospitals				<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(_numberOfRegions					)	<< "\t" << _numberOfRegions					<< std::endl;

	///// Parameters for runMCMCRGeneration function
	OUTPUT.MetaData << GET_VARIABLE_NAME(delta			)	<< "\t" << delta			<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(kRinit			)	<< "\t" << kRinit			<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(pDetectIndex	)	<< "\t" << pDetectIndex		<< std::endl;
	OUTPUT.MetaData << GET_VARIABLE_NAME(kIntro			)	<< "\t" << kIntro			<< std::endl;

	OUTPUT.MetaData.close(); 
}

std::string Choose_Simon_ScenarioName	(ModelRun &MR, double kIntro) //// unseen "arguments" are various global variables, e.g. _spatialLevel
{
	std::stringstream scenarioNameStream;
	scenarioNameStream << "V1.10correctTree";
	scenarioNameStream << "-data"		<< _dataSevereCases;
	scenarioNameStream << "-level"		<< _variabilityRAtClusterLevel;
	scenarioNameStream << "-spatial"	<< _spatialLevel;
	scenarioNameStream << "-heteroR"	<< _heterogeneityR;
	scenarioNameStream << "-reductionR" << _withReductionR;
	scenarioNameStream << "-seasonal"	<< _withSeasonality;
	scenarioNameStream << "-introOverD" << _introOverdispersed;
	if ((_introOverdispersed == 1) & (_estimOverdispersionIntro == 0)) scenarioNameStream << ".k" << kIntro;
	if (_estimOverdispersionIntro == 1) scenarioNameStream << ".estimKIntro";

	std::string Simon_scenarioName("");

	Simon_scenarioName = "";

	return Simon_scenarioName;
}
std::string Choose_InputFileName		(ModelRun &MR)
{
	string inputFileName;
	inputFileName = "MERS_forCpp.txt"; 
	return inputFileName; 
}

template <typename T>
std::string ToStringWithPrecision(const T a_value, const int n = 2)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}
std::string Choose_DJL_scenarioName	(ModelRun& MR)
{
	std::string DJL_ScenarioName = ""; 

	if (MR.Efficacy_Start == 0) 	DJL_ScenarioName = DJL_ScenarioName + "_NoHumans";  // i.e. if not vaccinating humans
	else							DJL_ScenarioName = DJL_ScenarioName + "_Eff_" + ToStringWithPrecision(MR.Efficacy_Start, 2); // i.e. if vaccinating humans
	
	if (MR.Efficacy_CamelControls > 0.0)
		DJL_ScenarioName = DJL_ScenarioName + "_CamelsControlEff_" + ToStringWithPrecision(MR.Efficacy_CamelControls, 2);
	
	if (MR.Efficacy_Start != 0) // i.e. if vaccinating humans
	{
		if (MR.VacCampStrategy == VaccCampaignStrategy::PROACTIVE)
		{
			DJL_ScenarioName = DJL_ScenarioName + "_ProAct";
			DJL_ScenarioName = DJL_ScenarioName + "_Dur_" + ToStringWithPrecision(MR.VaccineDuration		, 2);
			DJL_ScenarioName = DJL_ScenarioName + "_Lag_" + ToStringWithPrecision(MR.TimeSinceVaccination	, 2);
		}
		else if (MR.VacCampStrategy == VaccCampaignStrategy::REACTIVE)
		{
			if (MR.ReactLevel != ReactiveLevel::HOSPITAL)
			{
					 if (MR.ReactLevel == ReactiveLevel::REGIONAL) DJL_ScenarioName = DJL_ScenarioName + "_reg";
				else if (MR.ReactLevel == ReactiveLevel::NATIONAL) DJL_ScenarioName = DJL_ScenarioName + "_nat";
			}
		}

		if (MR.Coverage				!= 1.0	)	DJL_ScenarioName = DJL_ScenarioName + "_Cov"		+ ToStringWithPrecision		(MR.Coverage, 2			);
		if (MR.ImplementationDelay	!= 0	)	DJL_ScenarioName = DJL_ScenarioName + "_ImpDelay"	+ std::to_string			(MR.ImplementationDelay	);
		if (MR.ImmunityDelay		!= 0	)	DJL_ScenarioName = DJL_ScenarioName + "_VacDelay"	+ std::to_string			(MR.ImmunityDelay		);

			 if (MR.VaccinateAllHumans		)	DJL_ScenarioName = DJL_ScenarioName + "_Blanket"	;
		else if (MR.Vaccinate_HCW			)	DJL_ScenarioName = DJL_ScenarioName + "_vHCW"		;
	}

	if (!MR.RunOnCluster) DJL_ScenarioName = DJL_ScenarioName + "_Dummy";
	return DJL_ScenarioName; 
}

int Choose_numberOfHospitals(int dataSevereCases)
{
	int numberOfHospitals; 
		 if (_dataSevereCases == 0)	numberOfHospitals = 98;
	else if (_dataSevereCases == 1)	numberOfHospitals = 89;
	else							numberOfHospitals = 97;

	return numberOfHospitals; 
}

