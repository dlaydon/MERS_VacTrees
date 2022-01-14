
#include "Libraries.h"
#include "Macros.h"
#include "DataStructure.h"
#include "randlib_par.h"
#include "omp.h"
#include "DJL_Structs.h"
#include "Likelihood.h"
#include "Globals.h"
#include "InputOutput.h"
#include "MCMC.h"

long SEED;



void ChangeDownstreamCases(int CaseNumber, bool DeleteOrNot)
{
	//// function changes counterfactual _DeleteCase_CF container - DOES NOT CHANGE _secondaryCases container (as this would be actual, not counterfactual). 
	//// i.e. assign new value of DeleteOrNot to _DeleteCase_CF[CaseNumber]
	_DeleteCase_CF[CaseNumber] = DeleteOrNot;
	if (_secondaryCases[CaseNumber].empty()) return; //// if there are no secondary cases, can stop here... 
	//// ... but if not, loop through all secondary cases and amend their value of _DeleteCase_CF. And then recursively amend the values of _DeleteCase_CF for their secondary cases, and so on until no more secondary cases. 
	for (set<int>::iterator it = _secondaryCases[CaseNumber].begin(); it != _secondaryCases[CaseNumber].end(); ++it)
		ChangeDownstreamCases(*it, DeleteOrNot);
}
void FindAndDeleteCases()
{
	for (int CaseNumber = 0; CaseNumber < _numberOfCases; CaseNumber++)
		//if (_DeleteCase_CF[CaseNumber] == false) // only consider if statements and function calls below if CaseNumber not already deleted.
			if (_vaccinated[CaseNumber])
				if (_protected[CaseNumber])
					ChangeDownstreamCases(CaseNumber, /*DeleteOrNot =*/ true);
}



/////////////////		/////////////////		/////////////////		/////////////////		
/////////////////		MCMC / main functions
/////////////////		/////////////////		/////////////////		/////////////////		

bool DoYouVaccinate(ModelRun& MR, int &iInfectee, int Delay)
{
	return (MR.VacCampStrategy == VaccCampaignStrategy::PROACTIVE) ||
		(MR.ReactLevel == ReactiveLevel::HOSPITAL && _onset[iInfectee] >= _FirstOnsetInHosp		[_hospital[iInfectee]]	+ Delay) ||		//// if any delay would not have affected this patient. Note that because of MR.ImmunityDelay, people who could otherwise be vaccinated, but not mount immune response before they became a case, are considered unvaccinated. 
		(MR.ReactLevel == ReactiveLevel::REGIONAL && _onset[iInfectee] >= _FirstOnsetInRegion	[_region[iInfectee]]	+ Delay) ||
		(MR.ReactLevel == ReactiveLevel::NATIONAL && _onset[iInfectee] >= _FirstOnsetInCountry							+ Delay)	;
}

void runMCMC				(AllOutput &OUTPUT, ModelRun &MR, double optAR = 0.24, double delta = 0.)
{
	std::cout << "runMCMC: " << std::endl;
	fflush(stderr);	fflush(stdout);

	//// set column names 
	PutColNamesOn_AllOutput(OUTPUT);
	///// Initialize Output values. 
	InitializeOutputValues(OUTPUT);

	double currAR;
	int numberOfParameters = _parameter.size();
	vector<double>	numberOfMoveAccepted(numberOfParameters, 0.	);
	vector<int>		numberOfMoveProposed(numberOfParameters, 0	);
	double numberProposedUpdateInfector = 0., numberAcceptedUpdateInfector = 0.;

	double* _probaSource = new double[4];

	//// initialize _meanIndividualObsR
	for (int iInfectee = 1; iInfectee < _numberOfCases; iInfectee++)
		for (int H2Htype = 0; H2Htype < _numberOfH2Htype; H2Htype++)
			_meanIndividualObsR[iInfectee][H2Htype] = 0;

	//// initialize _logLik (global variable for log likelihood, not array of function pointers)
	_logLik = logLik();
	std::cout << setprecision(4) << "LL = " << _logLik << endl;

	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** 
	///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** ///// **** 
	///// **** ///// **** START MCMC
	
	std::cout << "================= START MCMC =================" << endl;
	fflush(stderr);	fflush(stdout);

	for (int iteration = 0; iteration < MR.NumIterations; ++iteration)
	{
		std::cout << iteration << ":";
		for (int thin = 0; thin < MR.StoreEvery; ++thin)
		{
			//// first update parameters (and adjust proposal step sizes based on current acceptance rates). 
			for (int parameter = 0; parameter < numberOfParameters; ++parameter)
				if (_rateForRandomWalk[parameter] > 0.)
				{
					numberOfMoveAccepted[parameter] += updateParameter(parameter);
					++numberOfMoveProposed[parameter];
					currAR = numberOfMoveAccepted[parameter] / numberOfMoveProposed[parameter];
					updateRate(_rateForRandomWalk[parameter], currAR, optAR, delta); //// updates _rateForRandomWalk[parameter]
				}

			//// Then update tree by adjusing infectors for all cases (infectees)
			for (int iInfectee = 1; iInfectee < _numberOfCases; ++iInfectee)
			{
				// Update infector
				numberAcceptedUpdateInfector += updateInfector(iInfectee); //// actual values of numberAcceptedUpdateInfector or numberProposedUpdateInfector are just diagnostics not used anywhere else in the program. updateInfector function is really a void function in that it changes LOTS of other global variables. 
				++numberProposedUpdateInfector;
			}

			///// updates part of likelihood 
			if (_variabilityRAtClusterLevel > vC_No_Hetero) updateExpR();

			//// updates w_n weights
			updateProbaOfTimingTrans();
		}

		//// Append Chains (logLik then parameters)
		if (OUTPUT.Write_Chains)
		{
			OUTPUT.Chains << _logLik << "\t";	
			for (int parameter = 0; parameter < numberOfParameters; ++parameter)	OUTPUT.Chains << _parameter[parameter] << "\t";
		}
		std::cout		<< _logLik << " ";

		//// Reset and recalculate _probaSource
		int H2Htype;
		for (int i = 0; i < 4; i++) _probaSource[i] = 0;
		for (int iInfectee = 1; iInfectee < _numberOfCases; iInfectee++)
			if (_infector[iInfectee] == AnimalReservoir) _probaSource[3]++;
			else 
			{
				H2Htype = computeH2Htype(iInfectee, _infector[iInfectee]);
				_probaSource[H2Htype]++;
			}
		for (int i = 0; i < 4; i++)
		{
			_probaSource[i] = _probaSource[i] / _numberOfCases;
			if (OUTPUT.Write_Chains) OUTPUT.Chains << _probaSource[i] << "\t";
		}

		int cumNbIntroducers = 0;
		for (int day = 0; day < _maxDay; day++) cumNbIntroducers += _numberOfIntroducers[day];
		if (OUTPUT.Write_Chains) OUTPUT.Chains << cumNbIntroducers << "\t";

		if (iteration > 0) simulateEpidemicWithCluster();
		else 
		{
			_simulNbIntro		= 0; _simulNbCase			= 0; _simulNbCluster			= 0; _simulMeanClusterSize			= 0; _simulMaxClusterSize = 0;
			_simulPClusterSize1 = 0; _simulPClusterSize10	= 0; _simulMaxNbCaseOverPeriod	= 0; _simulMaxNbClusterOverPeriod	= 0;
		}

		if (OUTPUT.Write_Chains)
		{
			OUTPUT.Chains << _simulNbIntro			<< "\t" << _simulNbCase			<< "\t" << _simulNbCluster				<< "\t" << _simulMeanClusterSize		<< "\t"		<< _simulMaxClusterSize << "\t";
			OUTPUT.Chains << _simulPClusterSize1	<< "\t" << _simulPClusterSize10 << "\t" << _simulMaxNbCaseOverPeriod	<< "\t" << _simulMaxNbClusterOverPeriod << std::endl;
		}
		
		
		//// **** //// **** //// **** //// **** //// **** //// **** /// **** //// **** //// **** //// **** //// **** //// **** /// **** //// **** //// **** //// **** //// **** //// **** 
		//// **** //// **** //// **** //// **** //// **** //// **** /// **** //// **** //// **** //// **** //// **** //// **** /// **** //// **** //// **** //// **** //// **** //// **** 
		/// **** //// **** //// **** //// **** //// **** //// **** COUNTERFACTUALS / TREE-PRUNING.

		/*	Make and output counterfactual tree

				i)		which cases have been vaccinated?
				ii)		which vaccinees are protected?
				iii)	which secondary cases wouldn't have happened? (in simplest case with no alternative infectors)
		*/

		for (int CFnum = 0; CFnum < MR.NumCFsPerTree; CFnum++)
		{
#ifdef PRINT_PROGRAM_PROGRESS
			std::cout << "reset, ";
			fflush(stdout); fflush(stderr);
#endif
			//// reset
			for (int iInfectee = 0; iInfectee < _numberOfCases; iInfectee++)
			{
				_protected		[iInfectee] = 0;
				_vaccinated		[iInfectee] = 0;
				_DeleteCase_CF	[iInfectee] = 0;
			}

#ifdef PRINT_PROGRAM_PROGRESS
			std::cout << "vaccinate, ";
			fflush(stdout); fflush(stderr); 
#endif

			// Camel control measures
			if (MR.Efficacy_CamelControls > 0.0)
				for (int iInfectee = 1; iInfectee < _numberOfCases; iInfectee++)
					if (_infector[iInfectee] == AnimalReservoir)
						if (BernoulliTrial(MR.Efficacy_CamelControls))
							ChangeDownstreamCases(iInfectee, /*DeleteOrNot =*/ true);

			//// vaccinate 
			if (MR.Efficacy_Start != 0)
			{
				if (MR.VaccinateAllHumans)
					for (int iInfectee = 1; iInfectee < _numberOfCases; iInfectee++)
						_vaccinated[iInfectee] = 1;
				else if (MR.Vaccinate_HCW)
				{
					for (int iInfectee = 1; iInfectee < _numberOfCases; iInfectee++)
						if (MR.Vaccinate_HCW && _HealthCareWorker[iInfectee])	//// vaccinate all healthcare workers (except reservoir)
							if (_HealthCareWorker[iInfectee])					//// if there was enough time for them to vaccinate 
								if (DoYouVaccinate(MR, iInfectee, MR.ImplementationDelay))
									_vaccinated[iInfectee] = BernoulliTrial(MR.Coverage);
				}
#ifdef PRINT_PROGRAM_PROGRESS
			std::cout << "protect, ";
			fflush(stdout); fflush(stderr);
#endif
				//// efficacy is:
				for (int iInfectee = 1; iInfectee < _numberOfCases; iInfectee++) 
					if (_vaccinated[iInfectee]) 
						if (DoYouVaccinate(MR, iInfectee, MR.ImplementationDelay + MR.ImmunityDelay)) //// if any delay would not have affected this patient. 
							_protected[iInfectee] = BernoulliTrial(MR.Efficacy_Current);

#ifdef PRINT_PROGRAM_PROGRESS
				std::cout << "delete cases, ";
				fflush(stdout); fflush(stderr);
#endif
			}

			//// Delete cases, and those downstream.
			FindAndDeleteCases();

			// // == == // // == == // // == == // // == == // // == == // // == == // // == == // // == == // // == == // // == == // // == == // // == == // // == == // // == == // // == == 
			// Record counterfactuals
			
			int NumCasesRemaining_CF		= _numberOfCases - 1;						// number of cases remaining counterfactual (minus 1 because of reservoir)
			int NumCases_HCW_Remaining_CF	= _numberOfHCWs;							// number of cases remaining in healthcare workers counterfactual
			int NumCases_nHCW_Remaining_CF	= NumCasesRemaining_CF - _numberOfHCWs;		// number of cases remaining in non-healthcare workers counterfactual
			int T_FinalCase_CF				= 0;										// counterfactual time of final case 
			int DeathsAverted				= 0;
			int DeathsAverted_HCW			= 0;
			int DeathsAverted_nHCW			= 0;

			for (int iInfectee = 1; iInfectee < _numberOfCases; iInfectee++)
			{
				if (_DeleteCase_CF[iInfectee])
				{
					NumCasesRemaining_CF--;
					if (_HealthCareWorker[iInfectee]) NumCases_HCW_Remaining_CF--; else NumCases_nHCW_Remaining_CF--;
					if (_Dead[iInfectee])
					{
						DeathsAverted++;
						if (_HealthCareWorker[iInfectee]) DeathsAverted_HCW++; else DeathsAverted_nHCW++;
					}
				}
				else
				{
					if (_onset[iInfectee] > T_FinalCase_CF) T_FinalCase_CF = _onset[iInfectee];

					//// add to epidemic curve. How many non-deleted cases were there on each day? For this iteration, add to this infectees onset date
					MR.CF_EpiCurves_Internal[iteration][_onset_week[iInfectee]]++;
					if (_Dead[iInfectee]) MR.CF_EpiCurves_Deaths_Internal[iteration][_onset_week[iInfectee]]++;
				}
			}

			MR.CF_Chains[iteration][MR.QIs.PropAverted			]	= ((double)_numberOfCases						- (double)NumCasesRemaining_CF		) / _numberOfCases	;	
			MR.CF_Chains[iteration][MR.QIs.PropAverted_HCW		]	= ((double)_numberOfHCWs						- (double)NumCases_HCW_Remaining_CF	) / _numberOfHCWs	;
			MR.CF_Chains[iteration][MR.QIs.PropAverted_nHCW		]	= ((double)(_numberOfCases - 1 - _numberOfHCWs)	- (double)NumCases_nHCW_Remaining_CF) / ((double)(_numberOfCases - 1 - _numberOfHCWs));
			MR.CF_Chains[iteration][MR.QIs.FinalCaseDate_CF		]	= T_FinalCase_CF;
			MR.CF_Chains[iteration][MR.QIs.DeathsAverted_CF		]	= DeathsAverted;
			MR.CF_Chains[iteration][MR.QIs.DeathsAverted_CF_HCW	]	= DeathsAverted_HCW;
			MR.CF_Chains[iteration][MR.QIs.DeathsAverted_CF_nHCW]	= DeathsAverted_nHCW;

#ifdef PRINT_PROGRAM_PROGRESS
			std::cout << "add to CR tree, ";
			fflush(stdout); fflush(stderr);
#endif 
			if (OUTPUT.Write_Trees_CF)	
				if (MR.OutputTreesEvery * int(iteration / MR.OutputTreesEvery) == iteration)
					for (int iInfectee = 1; iInfectee < _numberOfCases; iInfectee++)
						if (!_DeleteCase_CF[iInfectee]) //// i.e. if case still included after counterfactual campaign, then output everything that would normally be outputted. 
						{
							OUTPUT.Trees_CF << iteration << "\t" << CFnum << "\t" << iInfectee << "\t" << _hospital[iInfectee] << "\t" << _onset[iInfectee] << "\t";
							if (_infector[iInfectee] == AnimalReservoir) //// animal reservoir (coded zero throughout Cpp code, but -1 for R output). 
							{
								H2Htype = DifferentRegion;
								OUTPUT.Trees_CF << H2Htype << "\t" << -1 << "\t" << -1 << "\t" << -1 << endl;
							}
							else
							{
								H2Htype = computeH2Htype(iInfectee, _infector[iInfectee]);
								OUTPUT.Trees_CF << H2Htype << "\t" << _infector[iInfectee] << "\t" << _hospital[_infector[iInfectee]] << "\t" << _onset[_infector[iInfectee]] << endl;
							}
						}
		}



		if (MR.OutputTreesEvery * int(iteration / MR.OutputTreesEvery) == iteration)
		{
			if (OUTPUT.Write_SimulCluster)
				for (int iCluster = 0; iCluster < _simulNbCluster; iCluster++)
				{
					OUTPUT.SimulCluster << iteration << "\t" << _simulClusterStart[iCluster] << "\t" << _regionOfHospital[_simulClusterHospital[iCluster]] << "\t";
					OUTPUT.SimulCluster << _simulClusterHospital[iCluster] << "\t" << _simulNbInClusterBeforeDay[iCluster] << "\t" << _simulClusterR[iCluster][0] << endl;
				}

			//// Append Tree file (outputFileTree)
			if (OUTPUT.Write_Trees)	
				for (int iInfectee = 1; iInfectee < _numberOfCases; iInfectee++)
				{
					OUTPUT.Trees << iteration << "\t" << iInfectee << "\t" << _hospital[iInfectee] << "\t" << _region[iInfectee] << "\t" << _onset[iInfectee] << "\t" ;
					if (_infector[iInfectee] == AnimalReservoir) //// animal reservoir (coded zero throughout code, but -1 for R output). 
					{
						H2Htype = DifferentRegion;
						OUTPUT.Trees << H2Htype << "\t" << -1					<< "\t" << -1									<< "\t" << -1;
					} 
					else 
					{
						H2Htype = computeH2Htype(iInfectee, _infector[iInfectee]);
						OUTPUT.Trees << H2Htype << "\t" << _infector[iInfectee] << "\t" << _hospital[_infector[iInfectee]]	<< "\t" << _onset[_infector[iInfectee]]	;
					}
					OUTPUT.Trees << endl;
				}

		}

		if (OUTPUT.Write_ClusterR)
			if ((_variabilityRAtClusterLevel == vC_ClusterLevel) & (_heterogeneityR == 1))
			{
				for (int iCluster = 0; iCluster < _numberOfClusters; iCluster++) OUTPUT.ClusterR << _clusterExpR[iCluster][SameHosptial] << "\t";
				OUTPUT.ClusterR << endl;
			}
		if (OUTPUT.Write_ClusterReductionR)
			if (_heterogeneityReductionR == 1)
			{
				for (int iCluster = 0; iCluster<_numberOfClusters; iCluster++) OUTPUT.ClusterReductionR << _clusterReductionR[iCluster] << "\t";
				OUTPUT.ClusterReductionR << endl;
			}

		for (int iInfectee = 1; iInfectee < _numberOfCases; iInfectee++)
			for (int H2Htype = 0; H2Htype<_numberOfH2Htype; H2Htype++)    
				_meanIndividualObsR[iInfectee][H2Htype] += _individualObsR[iInfectee][H2Htype];

		std::cout << std::endl; 
		std::fflush(stderr); std::fflush(stdout);
	}
	//// MCMC finished

	std::cout << "================= END MCMC =================" << endl;
	fflush(stderr);	fflush(stdout);

	//// Process counterfactual epdemic curves
	if (OUTPUT.Write_CF_EpiCurves)
	{
		std::cout << "process CF epi curve" << std::endl; 
		double alpha				= 0.05;
		double RunningTotal			= 0;
		int NumPosteriorSamples		= MR.NumIterations - MR.BurnIn;
		int NumEntriesInTails		= floor	(alpha * NumPosteriorSamples / 2);
		int NumElementsForMedian	= NumPosteriorSamples / 2;
		int LowerIndex				= max(0, NumEntriesInTails - 1); 
		int UpperIndex				= min(NumPosteriorSamples - 1, NumPosteriorSamples - NumEntriesInTails - 1);
		std::vector<int> ThisDaysCFCasesOverAllIterations(NumPosteriorSamples, 0);

		std::cout << "NumIterations " << MR.NumIterations << " BurnIn " << MR.BurnIn << " NumPosteriorSamples " << NumPosteriorSamples << std::endl; 
		std::cout << "NumEntriesInTails " << NumEntriesInTails << " NumElementsForMedian " << NumElementsForMedian << " LowerIndex " << LowerIndex << " UpperIndex " << UpperIndex << std::endl;

		for (int Week = 0; Week < (MAX_ONSET_DAY / 7) + 1; Week++)
		{
			RunningTotal = 0;
			for (int iteration = MR.BurnIn; iteration < MR.NumIterations; iteration++) //// note start from BurnIn not zero
			{
				ThisDaysCFCasesOverAllIterations[iteration - MR.BurnIn] =  MR.CF_EpiCurves_Internal[iteration][Week];
				RunningTotal											+= MR.CF_EpiCurves_Internal[iteration][Week];
			}

			// sort this day's various iterations
			sort(ThisDaysCFCasesOverAllIterations.begin(), ThisDaysCFCasesOverAllIterations.end());

			// record
			MR.CF_EpiCurves[Week][MEAN]			= RunningTotal / NumPosteriorSamples;
			MR.CF_EpiCurves[Week][MEDIAN]		= ThisDaysCFCasesOverAllIterations[NumElementsForMedian];
			MR.CF_EpiCurves[Week][LOWER_CrI]	= ThisDaysCFCasesOverAllIterations[LowerIndex];
			MR.CF_EpiCurves[Week][UPPER_CrI]	= ThisDaysCFCasesOverAllIterations[UpperIndex];
		}
		std::cout << "Write CF_EpiCurves " << endl;

		// write column names
		OUTPUT.CF_EpiCurves << "Week" << "\t" << "Mean" << "\t" << "Median" << "\t" << "LowerCrI" << "\t" << "UpperCrI" << std::endl;
		for (int Week = 0; Week < (MAX_ONSET_DAY / 7) + 1; Week++)
			OUTPUT.CF_EpiCurves << Week << "\t" << MR.CF_EpiCurves[Week][MEAN] << "\t" << MR.CF_EpiCurves[Week][MEDIAN] << "\t" << MR.CF_EpiCurves[Week][LOWER_CrI] << "\t" << MR.CF_EpiCurves[Week][UPPER_CrI] << std::endl;
		std::cout << "CF_EpiCurves written " << endl;
	}

	//// Process counterfactual epdemic curves deaths 
	if (OUTPUT.Write_CF_EpiCurves_Deaths)
	{
		std::cout << "process CF epi curve - deaths" << std::endl; 
		double alpha				= 0.05;
		double RunningTotal			= 0;
		int NumPosteriorSamples		= MR.NumIterations - MR.BurnIn;
		int NumEntriesInTails		= floor	(alpha * NumPosteriorSamples / 2);
		int NumElementsForMedian	= NumPosteriorSamples / 2;
		int LowerIndex				= max(0, NumEntriesInTails - 1); 
		int UpperIndex				= min(NumPosteriorSamples - 1, NumPosteriorSamples - NumEntriesInTails - 1);
		std::vector<int> ThisDaysCFDeathsOverAllIterations(NumPosteriorSamples, 0);

		std::cout << "NumIterations " << MR.NumIterations << " BurnIn " << MR.BurnIn << " NumPosteriorSamples " << NumPosteriorSamples << std::endl; 
		std::cout << "NumEntriesInTails " << NumEntriesInTails << " NumElementsForMedian " << NumElementsForMedian << " LowerIndex " << LowerIndex << " UpperIndex " << UpperIndex << std::endl;

		for (int Week = 0; Week < (MAX_ONSET_DAY / 7) + 1; Week++)
		{
			RunningTotal = 0;
			for (int iteration = MR.BurnIn; iteration < MR.NumIterations; iteration++) //// note start from BurnIn not zero
			{
				ThisDaysCFDeathsOverAllIterations[iteration - MR.BurnIn]	=  MR.CF_EpiCurves_Deaths_Internal[iteration][Week];
				RunningTotal												+= MR.CF_EpiCurves_Deaths_Internal[iteration][Week];
			}

			// sort this day's various iterations
			sort(ThisDaysCFDeathsOverAllIterations.begin(), ThisDaysCFDeathsOverAllIterations.end());

			// record
			MR.CF_EpiCurves_Deaths[Week][MEAN]		= RunningTotal / NumPosteriorSamples;
			MR.CF_EpiCurves_Deaths[Week][MEDIAN]	= ThisDaysCFDeathsOverAllIterations[NumElementsForMedian];
			MR.CF_EpiCurves_Deaths[Week][LOWER_CrI]	= ThisDaysCFDeathsOverAllIterations[LowerIndex];
			MR.CF_EpiCurves_Deaths[Week][UPPER_CrI]	= ThisDaysCFDeathsOverAllIterations[UpperIndex];
		}
		std::cout << "Write CF_EpiCurves_Deaths " << endl;

		// write column names
		OUTPUT.CF_EpiCurves_Deaths << "Week" << "\t" << "Mean" << "\t" << "Median" << "\t" << "LowerCrI" << "\t" << "UpperCrI" << std::endl;
		for (int Week = 0; Week < (MAX_ONSET_DAY / 7) + 1; Week++)
			OUTPUT.CF_EpiCurves_Deaths << Week			<< "\t" << 
			MR.CF_EpiCurves_Deaths[Week][MEAN]			<< "\t" << 
			MR.CF_EpiCurves_Deaths[Week][MEDIAN]		<< "\t" << 
			MR.CF_EpiCurves_Deaths[Week][LOWER_CrI]		<< "\t" << 
			MR.CF_EpiCurves_Deaths[Week][UPPER_CrI]		<< std::endl;
		std::cout << "CF_EpiCurves_Deaths written " << endl;
	}



	//// write CF_Chains output
	if (OUTPUT.Write_CF_Chains)
	{
		std::cout << "Write CF_Chains " << endl; 
		for (int quant = 0; quant < MR.NumCFQuantities; quant++)
		{
			OUTPUT.CF_Chains << MR.CF_Names[quant]; 
			if (quant == MR.NumCFQuantities - 1) OUTPUT.CF_Chains << endl; else OUTPUT.CF_Chains << "\t";
		}

		for (int iteration = 0; iteration < MR.NumIterations; iteration++)
			for (int quant = 0; quant < MR.NumCFQuantities; quant++)
			{
				OUTPUT.CF_Chains << MR.CF_Chains[iteration][quant];
				if (quant == MR.NumCFQuantities - 1) OUTPUT.CF_Chains << endl; else OUTPUT.CF_Chains << "\t";
			}
		std::cout << "CF_Chains written " << endl;

	}
	std::fflush(stderr);	std::fflush(stdout);

	if (OUTPUT.Write_IndividualR) 
		for (int iInfectee = 1; iInfectee < _numberOfCases; iInfectee++)
		{
			OUTPUT.IndividualR << iInfectee << "\t" << _region[iInfectee] << "\t" << _hospital[iInfectee] << "\t" << _clusterOfCase[iInfectee] << "\t" << _onset[iInfectee] << "\t" << _rankOfCaseInCluster[iInfectee] << "\t";
			for (int H2Htype = 0; H2Htype<_numberOfH2Htype; H2Htype++)
			{
				_meanIndividualObsR[iInfectee][H2Htype] = _meanIndividualObsR[iInfectee][H2Htype] / MR.NumIterations;
				OUTPUT.IndividualR << _meanIndividualObsR[iInfectee][H2Htype] << "\t";
			}
			OUTPUT.IndividualR << endl;
		}

	cout << endl << "LL = " << _logLik << endl << "================== END MCMC ==================" << endl;
	if (OUTPUT.Write_Chains				) OUTPUT.Chains				.close();
	if (OUTPUT.Write_IndividualR		) OUTPUT.IndividualR		.close();
	if (OUTPUT.Write_ClusterR			) OUTPUT.ClusterR			.close();
	if (OUTPUT.Write_ClusterReductionR	) OUTPUT.ClusterReductionR	.close();
	if (OUTPUT.Write_SimulCluster		) OUTPUT.SimulCluster		.close();
	if (OUTPUT.Write_Trees				) OUTPUT.Trees				.close();
	if (OUTPUT.Write_Trees_CF			) OUTPUT.Trees_CF			.close();
	if (OUTPUT.Write_CF_Chains			) OUTPUT.CF_Chains			.close();
	if (OUTPUT.Write_CF_EpiCurves		) OUTPUT.CF_EpiCurves		.close();
	if (OUTPUT.Write_MetaData			) OUTPUT.MetaData			.close();

	std::cout << "runMCMC DONE" << std::endl;
	std::fflush(stderr);	std::fflush(stdout);
}

void runMCMCRGeneration		(AllOutput& OUTPUT, FileStrings_Struct FileStrings, double delta, double kRinit, double pDetectIndex, double kIntro, ModelRun &MR)
{
	std::cout << "runMCMCRGeneration " << std::endl;
	fflush(stderr);	fflush(stdout);

	///// Wrapper of loadAndInitializeData and runMCMC functions above, which specifies filepaths, initial parameter values, proposal values etc., then loads data before running MCMC

	// Parameters
	_beta = 0.05;

	if (_spatialLevel == SL_No_Space) _numberOfH2Htype = 1; //0: no space / 1: region / 2: cluster / 3: region+cluster
	else 
	{
		if (_spatialLevel < 3)	_numberOfH2Htype = 2;
		else					_numberOfH2Htype = 3;
	}

	//// **** //// **** //// **** Initial parameter values
	int numberOfParameters = 15;
	_parameter		= vector<double>(numberOfParameters, 0.);
	_parameter[0]	= 10;		// Mean generation time (taken from Lancet ID)
	_parameter[1]	= 4;		// Generation time std (taken from Lancet ID)
	_parameter[2]	= 1.0;		// R hosp
	_parameter[3]	= 1.0;		// R region
	_parameter[4]	= 1.0;		// R other region
	_parameter[5]	= kRinit;	// k hosp  //_cvParam=0: param5= k ; _cvParam=1: param5= cv=1/sqrt(k)
	_parameter[6]	= kRinit;	// k region
	_parameter[7]	= 0.1;		// Mean number of introductions at time 0
	_parameter[8]	= 0.;		// Mean number of introductions - exponential rate
	_parameter[9]	= kIntro;	// overdispersion intro
	_parameter[10]	= 0.7;		// reduction in cluster transmission over time
	_parameter[11]	= 0.5;		// power reduction in cluster transmission over time
	_parameter[12]	= 0.4;		// Time of year when seasonality peak (0=1 Jan), (1=31 Dec)
	_parameter[13]	= 0.5;
	_parameter[14]	= 1;


	//// **** //// **** //// **** Specify proposal step sizes / rates for random walk
	_rateForRandomWalk = vector<double>(numberOfParameters, 1.0);
	_rateForRandomWalk[0] = 0.1;																								// mean GT
	_rateForRandomWalk[1] = 0.1;																								// sd GT
	_rateForRandomWalk[2] = 0.1;																								// R level1
	_rateForRandomWalk[3] = (_spatialLevel > SL_No_Space) ? 0.1 : 0;															// R level2
	_rateForRandomWalk[4] = (_spatialLevel == SL_RegionAndHospital) ? 0.1 : 0;													// R level3
	_rateForRandomWalk[5] = ((_kRfix == 1) | (_variabilityRAtClusterLevel == vC_No_Hetero) | (_heterogeneityR == 0)) ? 0 : 0.1; // k level1
	_rateForRandomWalk[6] = ((_kRfix == 1) | (_heterogeneityRegionR == 0)) ? 0 : 0.1;											// k level1
	_rateForRandomWalk[7] = 1;																									// Mean number of introductions at time 0
	_rateForRandomWalk[8] = 0.005;																								// Mean number of introductions - exponential rate
	_rateForRandomWalk[9] = (_estimOverdispersionIntro == 0) ? 0 : 0.1;															// overdispersion
	_rateForRandomWalk[10] = (_withReductionR == 0) ? 0 : 0.1;																	//(_introOverdispersed == 0) ? 0:0.1; // overdispersion
	_rateForRandomWalk[11] = ((_withReductionR == 0) | (_heterogeneityReductionR == 0)) ? 0 : 0.1;								//(_introOverdispersed == 0) ? 0:0.1; // overdispersion
	_rateForRandomWalk[12] = (_withSeasonality>0) ? 0.1 : 0;
	_rateForRandomWalk[13] = (_withSeasonality>0) ? 0.1 : 0;
	_rateForRandomWalk[14] = 0;

	//// **** //// **** //// **** Specify parameter ranges / priors

	//// Lower limits
	_lowerLimit = vector<double>(numberOfParameters, -1e20);
	_lowerLimit[9] = -10.0;
	//_lowerLimit[13] = -1;

	//// Upper limits
	_upperLimit = vector<double>(numberOfParameters, 1e20);
	_upperLimit[0] = 30.0;
	_upperLimit[1] = 30.0;
	_upperLimit[2] = 30;
	_upperLimit[3] = 30;
	_upperLimit[4] = 30;
	_upperLimit[5] = (_cvParam == 0) ? 5 : 10000;
	_upperLimit[6] = (_cvParam == 0) ? 5 : 10000;
	_upperLimit[7] = 1e4;
	_upperLimit[8] = 10;
	_upperLimit[9] = 100;
	_upperLimit[10] = 10;
	_upperLimit[11] = 100;
	_upperLimit[12] = 1;
	_upperLimit[13] = 10;
	_upperLimit[14] = 1;

	_updatePointer		= vector<double(*)(int, double(*)())>(numberOfParameters, &updateParameterLogNormalProposal); //// _updatePointer is vector of pointers to functions. This line initializes all entries of this vector to equal updateParameterLogNormalProposal. Lines below change value of entries corresponding to parameters for which different updating functions are used. 
	_updatePointer[0]	= &updateParameterGT;
	_updatePointer[1]	= &updateParameterGT;
	_updatePointer[8]	= &updateParameterNormalProposal;

	_logLikPointer = vector<double(*)()>(numberOfParameters, &logLik); //// Logic similar to _updatePointer above. _logLikPointer vector of functions to calculate log likelihood. All entries are initialized to logLik

	// Load data
	loadAndInitializeData(FileStrings.inputFile, MR, false);

	// Run (and time) MCMC simulation
	clock_t tBeg = clock();
	runMCMC(OUTPUT, MR, 0.24, delta);
	clock_t tEnd = clock();
	double elapsed_secs = double((double)tEnd - (double)tBeg) / CLOCKS_PER_SEC;
	cout << endl << "This simulation took: " << elapsed_secs << "s" << endl << endl;
	cout << "--------- Parameters and MCMC Rates ----------------" << endl;
	for (int parameter = 0; parameter < numberOfParameters; ++parameter)
		cout << parameter << ": " << fixed << setw(15) << setprecision(4) << _parameter[parameter] << ", " << _rateForRandomWalk[parameter] << endl;
	cout << "----------------------------------------------------" << endl;
	std::fflush(stderr);	std::fflush(stdout);
}

void CorrectAndProcessInputParams(ModelRun& MR)
{
	std::cout << "CorrectAndProcessInputParams " << std::endl;

	if (MR.VacCampStrategy == VaccCampaignStrategy::REACTIVE)
	{
		MR.Efficacy_Current		= MR.Efficacy_Start;
		MR.TimeSinceVaccination = 0.0;
	}
	else if (MR.VacCampStrategy == VaccCampaignStrategy::PROACTIVE)
	{
		if (MR.VaccineDuration == 0.0) //// i.e. infinite duration
			MR.Efficacy_Current = MR.Efficacy_Start; 
		else
			MR.Efficacy_Current = MR.Efficacy_Start * exp(-MR.TimeSinceVaccination / MR.VaccineDuration); // set Efficacy_Current (can extend to have person specific TimeSinceVaccination if required). 
	}
	else std::cout  << "CorrectAndProcessInputParams error: MR.VacCampStrategy not recognized." << std::endl;
}

int main(int argc, char *argv[])
{
	//// wrapper of runMCMCRGeneration (which itself is a wrapper of runMCMC function), that sets up various "scenario" variables 
	ModelRun MR; 
	if (MR.RunOnCluster == false)
		std::cout << "MR.RunOnCluster == false" << std::endl;
	fflush(stderr);	fflush(stdout);

	if (MR.RunOnCluster)
	{
		string pParamFileName = argv[1];
		ReadInParams(MR, pParamFileName);
	}
	fflush(stderr);	fflush(stdout);
	MR.init();
	CorrectAndProcessInputParams(MR); 
	fflush(stderr);	fflush(stdout);

	std::cout << "Strategy: " << Convert_VaccCampaignStrategy_FromEnumClass(MR.VacCampStrategy) << std::endl; 
	std::cout << "Efficacy_Start " << MR.Efficacy_Start << ", Efficacy_Current " << MR.Efficacy_Current << std::endl;
	if (MR.VacCampStrategy == VaccCampaignStrategy::PROACTIVE)
		std::cout << "VaccineDuration " << MR.VaccineDuration << ", TimeSinceVaccination " << MR.TimeSinceVaccination << std::endl;

	//// Examples of scenario variables include _withReductionR, _variabilityRAtClusterLevel etc. 
	//// Input //// //// //// //// //// 
	_dataSevereCases	= 0;				// 0: all cases			; 1: severe cases										; 2: all cases no cluster 44
	_withReductionR		= 1;				// 0: no reduction R	; 1: reduction with cumulated number of clusters		; 2: reduction with delay since start of cluster
	_withSeasonality	= No_Seasonality;	// 0: no seasonality	; 1: seasonality (effect seasonality can't be below 1)	; 2: true cos function
	_introOverdispersed = 0;				// 0: no				; 1: yes
	
	//// Initialize seeds
	initSeeds(MR.seed1, MR.seed2, MR.max_threads);  // counterfactuals
	SEED = -194837829;								// MCMC

	//-----------
	double kIntro = 0.2;
	_estimOverdispersionIntro	= 0;
	_variabilityRAtClusterLevel = vC_ClusterLevel;			// 0: no			/ 1: individual		/ 2: cluster
	_spatialLevel				= SL_RegionAndHospital;		// 0: no space		/ 1: region			/ 2: hospital	/ 3: region+hospital
	_heterogeneityR				= 1;
	_heterogeneityReductionR	= 0;

	/////////////////
	_maxDuration			= 21;
	_trackGeneration		= 0;
	double pDetectIndex		= 0;
	_kRfix					= 0;
	_heterogeneityRegionR	= 0;
	_cvParam				= 1;
	double kRinit			= 1;
	if (_cvParam == 1) kRinit = 1. / pow(kRinit, 0.5);
	double delta			= 0;

	//// Choose ScenarioName and input file name. 
	string Simon_scenarioName		= Choose_Simon_ScenarioName	(MR, kIntro)	; 	std::cout << "Simon_scenarioName = "	<< Simon_scenarioName	<< std::endl;
	string inputFileName			= Choose_InputFileName		(MR)			; 	std::cout << "inputFileName = "			<< inputFileName		<< std::endl;
	string DJLScenarioName			= Choose_DJL_scenarioName	(MR)			; 	std::cout << "DJLScenarioName = "		<< DJLScenarioName		<< std::endl;

	//// AssignNumber of Hospitals & regions (global variables)
	_numberOfHospitals	= Choose_numberOfHospitals(_dataSevereCases); 
	_numberOfRegions	= 11;

	//// File paths
	FileStrings_Struct FileStrings; 
	FileStrings.init (MR.RunOnCluster, Simon_scenarioName, inputFileName, DJLScenarioName);
	AllOutput OUTPUT;
	OUTPUT.init(MR.RunOnCluster, FileStrings);

	///// Write model meta data. 

	if (OUTPUT.Write_MetaData)	WriteModelMetaData(OUTPUT, MR, FileStrings, delta, kRinit, pDetectIndex, kIntro);
	runMCMCRGeneration(OUTPUT, FileStrings, delta, kRinit, pDetectIndex, kIntro, MR);

	return 0;
}
