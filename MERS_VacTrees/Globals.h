#pragma once

// data
extern int					_numberOfCases; // Number of cases + community (case 0)
extern int					_numberOfCases_PruningSubset;	// Number of cases + community (case 0) between MR.MinDay_Pruning and MR.MaxDay_Pruning
extern int					_numberOfDeaths_PruningSubset;	// Number of cases + community (case 0) between MR.MinDay_Pruning and MR.MaxDay_Pruning
extern int					_numberOfHospitals;
extern int					_numberOfRegions;
extern vector<int>			_regionOfHospital;
extern vector<int>			_numberOfHospitalsPerRegion;
extern vector<vector<int>>  _hospitalInRegion;
extern int					_numberOfHCWs; // number of healthcare workers

extern int _minDay, _maxDay;
extern vector<int>				_nameCase;
extern vector<int>				_onset;
extern vector<int>				_onset_week;
extern vector<int>				_hospital;
extern vector<int>				_region;
extern vector<int>				_infector;				// indexed by i) case; _infector[infectee] gives the infector of infectee
extern int						_numberOfH2Htype;		// 0: hospital; 1: region; 2: other region
extern vector<vector<int>>		_individualObsR;		/// Observed R0 by case and definition IN PARTICULAR PARAMATER SAMPLE,  indexed by i) case; ii) H2H type
extern vector<vector<double>>	_meanIndividualObsR;	/// Observed R0 by case and definition in MEAN. indexed by i) case; ii) H2H type
extern vector<vector<double>>	_individualExpR;
extern vector<vector<double>>	_clusterExpR;			//// within cluster EXPECTED reproduction numbers (drawn from Gamma dist). indexed by i) cluster; ii) h2h type (although second index always set to zero = SameHospital so this is a relic from older analyses)
extern vector<double>			_clusterReductionR;
extern vector<double>		_regionExpR;
extern int					_heterogeneityRegionR;

extern map<int, int>		_numberOfIntroducers;			// [day] Number of introducer on that day
extern int					_numberOfIntroducersWith0Sec;	// number of introducers with zero secondary infectees
extern vector<int>			_generation;					// [iCase] generation of case iCase. Think this means "how many infections were there between infectee and animal reservoir?" If one then 0th generation, if two then 1st generation ... 
extern vector<set<int>>		_secondaryCases;				// _secondaryCases[infector] = set of infected by infector
extern vector<set<int>>		_possibleInfectors;				// set in loadAndInitializeData function, and then never altered. 

extern int					_numberOfClusters;
extern int					_minDelayBetweenClusters;
extern vector<int>			_clusterOfCase;
extern vector<int>			_rankOfCaseInCluster;			// Num cases in cluster prior to onset of this case
extern vector<int>			_startOfCluster;
extern vector<int>			_hospitalOfCluster;
extern vector<set<int>>		_casesInCluster;			//// set of all cases in particular cluster.	indexed by i) cluster.	Set in loadAndInitializeData and never changed. 
extern vector<set<int>>		_casesInRegion;				//// set of all cases in particular region.		indexed by i) region.	Set in loadAndInitializeData and never changed. 

//simul
extern vector <int>			_simulCaseID;
extern vector <int>			_simulOnset;
extern vector <int>			_simulInfectorID;
extern vector <int>			_simulHospID;
extern vector <int>			_simulRegionID;
extern vector <vector<int>>	_simulNbSecCases;
extern vector<set<int>>		_simulCaseOnDay;

extern vector <int>			_simulClusterOfCase;
extern vector <vector<double>> _simulClusterR;
extern vector <int>			_simulClusterStart;
extern vector <int>			_simulClusterHospital;
extern vector <int>			_simulNbInClusterBeforeDay;
extern vector <int>			_simulNbInClusterOnDay;
extern vector <int>			_simulLastOnsetInHospital;
extern vector <int>			_simulLastClusterInHospital;
extern vector <double>			_simulRegionR;

extern int		_simulNbIntro;
extern int		_simulNbCase;
extern int		_simulNbCluster;
extern double	_simulMeanClusterSize;
extern double	_simulPClusterSize1;
extern double	_simulPClusterSize10;
extern int		_simulMaxClusterSize;
extern int		_simulMaxNbCaseOverPeriod;
extern int		_simulMaxNbClusterOverPeriod;

//parameters
extern vector<double> _parameter;
extern vector<double> _rateForRandomWalk;
extern vector<double> _lowerLimit;
extern vector<double> _upperLimit;
extern vector<double (*)(int, double (*)())> _updatePointer;
extern vector<double (*)()> _logLikPointer;
extern double _logLik;

extern int		_maxDuration;
extern double	_beta; // Community
extern int		_trackGeneration;

extern vector<double>			_densityGT; // Generation time density
extern vector<vector<double>>	_probaOfTimingTrans;				//// un-normalised w_n weights for each infectee/infector (i.e. from page 2 of supporting info (section "Update of the source of infection of case i"))
extern vector<vector<double>>	_weightedProbaOfTimingTrans;		////    normalised w_n weights for each infectee/infector (i.e. from page 2 of supporting info (section "Update of the source of infection of case i"))

extern int _variabilityRAtClusterLevel;	// 0: no heterogeneity / 1: individual level /2: cluster level
extern int _heterogeneityR;
extern int _heterogeneityReductionR;
extern int _withReductionR;				// 0: no reduction; 1: with cum nb cases; 2: with delay since start cluster
extern int _withSeasonality;			// 0: no seasonality	; 1: seasonality (effect seasonality can't be below 1)	; 2: true cos function
extern int _dataSevereCases;


extern int _spatialLevel;					// 0: no space / 1:region / 2: hospital/ 3: region+hospital
extern int _introOverdispersed;				// 0: no overdispersion;  1: overdispersion
extern int _kRfix;
extern int _cvParam;						// 0: param5 = k ; 1: param5 = cv = 1/sqrt(k)
extern int _estimOverdispersionIntro;

//// Containers (anything with suffix _CF stands for counterfactual). 
// in simplest case, will assume 100% efficacy, so that _vaccinated and _protected are the same. 
extern vector<bool>			_vaccinated; 			// indexed by i) case;
extern vector<bool>			_protected; 			// indexed by i) case;
extern vector<bool>			_withinPruningWindow;	// indexed by i) case; Was subject's onset date within pruning window?
extern vector<bool>			_DeleteCase_CF;			// indexed by i) case; 
extern vector<int>			_HealthCareWorker;		// Healthcare worker. indexed by i) case; 0 = Not HCW			; 1 = HCW
extern vector<int>			_Dead;					// indexed by i) case; 					  0 = Alive				; 1 = Dead
extern vector<int>			_FirstOnsetInHosp;		// time of first onset at hospital. Indexed by i) hospital. 
extern vector<int>			_FirstOnsetInRegion; 	// time of first onset at region. Indexed by i) region. 
extern int					_FirstOnsetInCountry; 	// time of first onset in country. 
extern vector<set<int>>		_secondaryCases_CF;		// Counterfactual: _secondaryCases[infector] = set of infected by infector. 
extern vector<double>		_EfficacyCurrent;		// indexed by i) case. What is the remaining/residual efficacy for this case at the time of their symptom onset.
extern vector<int>			_DayTriggerReached;		// Date that trigger reached. initialized to MAX_ONSET_DAY + 1 (i.e. not reached). Indexed by i) either hopsital or region (or simply country) depending on whether reacting at hospital, regional or national level.
extern vector<vector<int>>	_EpiCurves;				// indexed by i) hospital or region; ii) day. Used only when doing triggers.
