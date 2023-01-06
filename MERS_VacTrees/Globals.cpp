#include "Libraries.h"

// data
int					_numberOfCases; // Number of cases + community (case 0)
int					_numberOfCases_PruningSubset; // Number of cases + community (case 0) between MR.MinDay_Pruning and MR.MaxDay_Pruning
int					_numberOfDeaths_PruningSubset; // Number of cases + community (case 0) between MR.MinDay_Pruning and MR.MaxDay_Pruning
int					_numberOfHospitals;
int					_numberOfRegions;
vector<int>			_regionOfHospital;
vector<int>			_numberOfHospitalsPerRegion;
vector<vector<int>> _hospitalInRegion;
int					_numberOfHCWs; // number of healthcare workers

int _minDay, _maxDay;
vector<int>				_nameCase;
vector<int>				_onset;
vector<int>				_onset_week;
vector<int>				_hospital;
vector<int>				_region;
vector<int>				_infector;				// indexed by i) case; _infector[infectee] gives the infector of infectee
int						_numberOfH2Htype;		// 0: hospital; 1: region; 2: other region
vector<vector<int>>		_individualObsR;		/// Observed R0 by case and definition IN PARTICULAR PARAMATER SAMPLE,  indexed by i) case; ii) H2H type
vector<vector<double>>	_meanIndividualObsR;	/// Observed R0 by case and definition in MEAN. indexed by i) case; ii) H2H type

vector<vector<double>>	_individualExpR;
vector<vector<double>>	_clusterExpR;			//// within cluster EXPECTED reproduction numbers (drawn from Gamma dist). indexed by i) cluster; ii) h2h type (although second index always set to zero = SameHospital - relic from older analyses)
vector<double>			_clusterReductionR;

vector<double>		_regionExpR;
int					_heterogeneityRegionR;

map<int, int>		_numberOfIntroducers;			// [day] Number of introducer on that day
int					_numberOfIntroducersWith0Sec;	// Number of introducers with zero secondary infectees
vector<int>			_generation;					// [iCase] generation of case iCase, i.e. "how many infections were there between infectee and animal reservoir?" If one then 0th generation, if two then 1st generation ... 
vector<set<int>>	_secondaryCases;				// _secondaryCases[infector] = set of infected by infector
vector<set<int>>	_possibleInfectors;				// set in loadAndInitializeData function, and then never altered. 

int					_numberOfClusters;
int					_minDelayBetweenClusters;
vector<int>			_clusterOfCase;
vector<int>			_rankOfCaseInCluster;		//// Num cases in cluster prior to onset of this case
vector<int>			_startOfCluster;
vector<int>			_hospitalOfCluster;
vector<set<int>>	_casesInCluster;			//// set of all cases in particular cluster.	indexed by i) cluster.	Set in loadAndInitializeData and never changed. 
vector<set<int>>	_casesInRegion;				//// set of all cases in particular region.		indexed by i) region.	Set in loadAndInitializeData and never changed. 

vector <int>			_simulCaseID;
vector <int>			_simulOnset;
vector <int>			_simulInfectorID;
vector <int>			_simulHospID;
vector <int>			_simulRegionID;
vector <vector<int>>	_simulNbSecCases;
vector<set<int>>		_simulCaseOnDay;

vector <int>			_simulClusterOfCase;
vector <vector<double>> _simulClusterR;
vector <int>			_simulClusterStart;
vector <int>			_simulClusterHospital;
vector <int>			_simulNbInClusterBeforeDay;
vector <int>			_simulNbInClusterOnDay;
vector <int>			_simulLastOnsetInHospital;
vector <int>			_simulLastClusterInHospital;
vector <double>			_simulRegionR;

int		_simulNbIntro;
int		_simulNbCase;
int		_simulNbCluster;
double	_simulMeanClusterSize;
double	_simulPClusterSize1;
double	_simulPClusterSize10;
int		_simulMaxClusterSize;
int		_simulMaxNbCaseOverPeriod;
int		_simulMaxNbClusterOverPeriod;

//parameters
vector<double> _parameter;
vector<double> _rateForRandomWalk;
vector<double> _lowerLimit;
vector<double> _upperLimit;
vector<double (*)(int, double (*)())> _updatePointer;
vector<double (*)()> _logLikPointer;
double _logLik;

int		_maxDuration;
double	_beta; // Community
int		_trackGeneration;

vector<double>			_densityGT;							// Generation time density
vector<vector<double>>	_probaOfTimingTrans;				//// un-normalised w_n weights for each infectee/infector (i.e. from page 2 of PNAS supporting info (section "Update of the source of infection of case i"))
vector<vector<double>>	_weightedProbaOfTimingTrans;		////    normalised w_n weights for each infectee/infector (i.e. from page 2 of PNAS supporting info (section "Update of the source of infection of case i"))

int _variabilityRAtClusterLevel;	// 0: no heterogeneity / 1: individual level /2: cluster level
int _heterogeneityR;
int _heterogeneityReductionR;
int _withReductionR;				// 0: no reduction; 1: with cum nb cases; 2: with delay since start cluster
int _withSeasonality;				// 0: no seasonality	; 1: seasonality (effect seasonality can't be below 1)	; 2: true cos function
int _dataSevereCases;

int _spatialLevel;					// 0: no space / 1:region / 2: hospital/ 3: region+hospital
int _introOverdispersed;			// 0: no overdispersion;  1: overdispersion
int _kRfix;
int _cvParam;						// 0: param5 = k ; 1: param5 = cv = 1/sqrt(k)
int _estimOverdispersionIntro;

//// Containers (anything with suffix _CF stands for counterfactual). 
vector<bool>		_vaccinated; 			// indexed by i) case;
vector<bool>		_withinPruningWindow;	// indexed by i) case; Was subject's onset date within pruning window?
vector<bool>		_protected; 			// indexed by i) case;
vector<bool>		_DeleteCase_CF;			// indexed by i) case; 
vector<int>			_HealthCareWorker;		// Healthcare worker. indexed by i) case; 0 = Not HCW			; 1 = HCW
vector<int>			_Dead;					// indexed by i) case; 					  0 = Alive				; 1 = Dead
vector<int>			_FirstOnsetInHosp;		// time of first onset at hospital. Indexed by i) hospital. 
vector<int>			_FirstOnsetInRegion; 	// time of first onset at reion. Indexed by i) region. 
int					_FirstOnsetInCountry; 	// time of first onset in country. 
vector<set<int>>	_secondaryCases_CF;		// Counterfactual: _secondaryCases[infector] = set of infected by infector. 
vector<double>		_EfficacyCurrent;		// indexed by i) case. What is the remaining/residual efficacy for this case at the time of their symptom onset.
vector<int>			_DayTriggerReached;		// Date that trigger reached. initialized to MAX_ONSET_DAY + 1 (i.e. not reached). Indexed by i) either hopsital or region (or simply country) depending on whether reacting at hospital, regional or national level.
vector<vector<int>>	_EpiCurves;				// indexed by i) hospital or region; ii) day. Used only when doing triggers.
