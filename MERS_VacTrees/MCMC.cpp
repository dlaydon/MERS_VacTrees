#include "Libraries.h"
#include "Macros.h"
#include "DataStructure.h"
#include "randlib_par.h"
#include "omp.h"
#include "Structs.h"
#include "Likelihood.h"
#include "MCMC.h"
#include "Globals.h"
#include "InputOutput.h"


/////////////////		/////////////////		/////////////////		/////////////////		
/////////////////		Parameter update functions
/////////////////		/////////////////		/////////////////		/////////////////		

//// The doubles returned here aren't very important. Better to think of these as voids, in that they change important stuff like logLik, _parameter values etc.
//// Also, all parameters updated using updateParameter function / _updatePointer vector of functions. Each parameter updated using 1 of 3 options: i) updateParameterNormalProposal; ii) updateParameterLogNormalProposal; iii) updateParameterGT
double updateParameterNormalProposal	(int parameterNumber, double (*logLik)())
{
	double oldValue		= _parameter[parameterNumber];
	double newValue		= oldValue + _rateForRandomWalk[parameterNumber] * rnorm();
	double logProposal	= 0.;	//// log of the correction factor

	//// reject if proposed value outside parameter ranges
		 if (newValue < _lowerLimit[parameterNumber]) return 0.;
	else if (newValue > _upperLimit[parameterNumber]) return 0.;

	///// Change relevant parameter of _parameter vector
	_parameter[parameterNumber] = newValue;
	double newLogLik			= logLik(); //// recalculate logLikelihood value in full
	double Q					= newLogLik - _logLik + logProposal;

	if (log(runif()) < Q)	///// i.e. if you accept
	{
		_logLik = newLogLik;	//// reset _logLik
		return 1.;				//// return 1 which			ADDS		to running total of number of acceptances when this function is called 
	}
	else 					///// i.e. if you reject
	{
		_parameter[parameterNumber] = oldValue;
		return 0.;				//// return 0 which		DOES NOT ADD	to running total of number of acceptances when this function is called 
	}
}
double updateParameterLogNormalProposal	(int parameterNumber, double (*logLik)())
{
	double oldValue		= _parameter[parameterNumber];
	double exponent		= _rateForRandomWalk[parameterNumber] * rnorm();
	double newValue		= oldValue * exp(exponent);
	double logProposal	= exponent; // log(newValue) - log(oldValue); 	//// log of the correction factor

	//// reject if proposed value outside parameter ranges
		 if (newValue < _lowerLimit[parameterNumber]) return 0.;
	else if (newValue > _upperLimit[parameterNumber]) return 0.;

	///// Change relevant parameter of _parameter vector
	_parameter[parameterNumber] = newValue;
	double newLogLik			= logLik(); //// recalculate logLikelihood value in full 
	double Q					= newLogLik - _logLik + logProposal;

	if (log(runif()) < Q)	///// i.e. if you accept
	{
		_logLik = newLogLik;	//// reset _logLik
		return 1.;				//// return 1 which			ADDS		to running total of number of acceptances when this function is called 
	}
	else 					///// i.e. if you reject
	{
		_parameter[parameterNumber] = oldValue;
		return 0.;				//// return 0 which		DOES NOT ADD	to running total of number of acceptances when this function is called 
	}
}
double updateParameterGT				(int parameterNumber, double(*logLik)())
{
	double oldValue		= _parameter[parameterNumber];
	double exponent		= _rateForRandomWalk[parameterNumber] * rnorm();
	double newValue		= oldValue * exp(exponent);
	double logProposal	= exponent; // log(newValue) - log(oldValue);

	//// reject if proposed value outside parameter ranges
	if (newValue > _upperLimit[parameterNumber]) return 0.;

	///// Change relevant parameter of _parameter vector 
	_parameter[parameterNumber]	= newValue;
	updateDensityGT();				//// updates _densityGT vector (probability density for generation time).

	double newLogLik			= logLik(); //// recalculate logLikelihood value in full
	double Q					= newLogLik - _logLik + logProposal;

	if (log(runif()) < Q)	///// i.e. if you accept
	{
		_logLik = newLogLik;	//// reset _logLik
		return 1.;				//// return 1 which			ADDS		to running total of number of acceptances when this function is called 
	}
	else  					///// i.e. if you reject
	{
		_parameter[parameterNumber] = oldValue;
		updateDensityGT();		//// updates _densityGT vector (probability density for generation time)
		return 0.;				//// return 0 which		DOES NOT ADD	to running total of number of acceptances when this function is called 
	}
}

//// Depending on the parameter number, this calls either of the three functions above. 
double updateParameter					(int parameterNumber) 
{
	return _updatePointer[parameterNumber](parameterNumber, _logLikPointer[parameterNumber]);
}

//// updateExpR not an paramter update function per se, as it assumes parameters have already gone through an ACCEPT/REJECT step.
///// Instead, updateExpR updates varrious quantities associated with the likelihood and that are necessary to calculate the likelihood, e.g. the expected R's at cluster level (or individual level depending on _variabilityRAtClusterLevel)
double updateExpR						()	
{
	//// now that parameters and infectors have been updated, need to update various quantities and total likelihood. 
	int H2Htype = SameHosptial;
	double oldLL, newLL;

	if (_variabilityRAtClusterLevel == vC_IndividualLevel)
		for (int iInfector = 1; iInfector < _numberOfCases; ++iInfector)
		{
			oldLL = ldgamma(_individualExpR[iInfector][H2Htype], aR(H2Htype), bR(H2Htype));
			oldLL += ldensityIndividualR(iInfector, H2Htype);

			_individualExpR[iInfector][H2Htype] = rgamma(aR(H2Htype) + _individualObsR[iInfector][H2Htype], bR(H2Htype) + 1);

			newLL = ldgamma(_individualExpR[iInfector][H2Htype], aR(H2Htype), bR(H2Htype));
			newLL += ldensityIndividualR(iInfector, H2Htype);

			_logLik += newLL - oldLL;
		}
	else if (_variabilityRAtClusterLevel == vC_ClusterLevel)
	{
		double sumObsRCluster, sumReductionRCluster, nbCasesCluster;
		if (_heterogeneityR == 1)
			for (int iCluster = 0; iCluster < _numberOfClusters; ++iCluster)
			{
				//// identify and calculate contribution to total likelhood from this cluster, BEFORE expected R for this cluster (_clusterExpR) changed. 
				oldLL					= ldgamma(_clusterExpR[iCluster][H2Htype], aR(H2Htype), bR(H2Htype)); //// L^cluster
				sumObsRCluster			= 0;
				nbCasesCluster			= 0;
				sumReductionRCluster	= 0;
				for (set<int>::iterator it = _casesInCluster[iCluster].begin(); it != _casesInCluster[iCluster].end(); ++it)
				{
					nbCasesCluster++;
					sumObsRCluster			+= _individualObsR[*it][H2Htype];
					sumReductionRCluster	+= reductionRForCase((*it));
					oldLL					+= ldensityIndividualR(*it, H2Htype); //// first line of L_n^trans definition (summed / multiplied over case numbers n in this cluster). 
				}

				_clusterExpR[iCluster][H2Htype] = rgamma(aR(H2Htype) + sumObsRCluster, bR(H2Htype) + sumReductionRCluster);

				//// identify and calculate contribution to total likelhood from this cluster, AFTER expected R for this cluster (_clusterExpR) changed. 
				newLL = ldgamma(_clusterExpR[iCluster][H2Htype], aR(H2Htype), bR(H2Htype));
				for (set<int>::iterator it = _casesInCluster[iCluster].begin(); it != _casesInCluster[iCluster].end(); ++it)
					newLL += ldensityIndividualR(*it, H2Htype); //// first line of L_n^trans definition (summed / multiplied over case numbers n in this cluster). 

				_logLik += newLL - oldLL;
			}

		if (_heterogeneityReductionR == 1)
			for (int iCluster = 0; iCluster < _numberOfClusters; ++iCluster)
			{
				oldLL = ldgamma(_clusterReductionR[iCluster], aReductionR(), bReductionR());
				for (set<int>::iterator it = _casesInCluster[iCluster].begin(); it != _casesInCluster[iCluster].end(); ++it)
					oldLL += ldensityIndividualR(*it, H2Htype);

				double oldValue = _clusterReductionR[iCluster];
				double exponent = 1 * rnorm();
				double newValue = oldValue * exp(exponent);
				double logProposal = exponent; // log(newValue) - log(oldValue);

				_clusterReductionR[iCluster] = newValue;

				newLL = ldgamma(_clusterReductionR[iCluster], aReductionR(), bReductionR());
				for (set<int>::iterator it = _casesInCluster[iCluster].begin(); it != _casesInCluster[iCluster].end(); ++it)
					newLL += ldensityIndividualR(*it, H2Htype);

				double Q = newLL - oldLL + logProposal;
				if (log(runif()) < Q)	_logLik += newLL - oldLL;
				else					_clusterReductionR[iCluster] = oldValue;
			}
	}

	H2Htype = SameRegion;
	int iRegion;
	if (_heterogeneityRegionR == 1)
		for (iRegion = 0; iRegion < _numberOfRegions; iRegion++)
		{
			oldLL = ldgamma(_regionExpR[iRegion], aR(H2Htype), bR(H2Htype));
			for (set<int>::iterator it = _casesInRegion[iRegion].begin(); it != _casesInRegion[iRegion].end(); ++it)
				oldLL += ldensityIndividualR(*it, H2Htype);

			double oldValue		= _regionExpR[iRegion];
			double exponent		= 1 * rnorm();
			double newValue		= oldValue * exp(exponent);
			double logProposal	= exponent; // log(newValue) - log(oldValue);

			_regionExpR[iRegion] = newValue;

			newLL = ldgamma(_regionExpR[iRegion], aR(H2Htype), bR(H2Htype));
			for (set<int>::iterator it = _casesInRegion[iRegion].begin(); it != _casesInRegion[iRegion].end(); ++it)
				newLL += ldensityIndividualR(*it, H2Htype);

			double Q = newLL - oldLL + logProposal;
			if (log(runif()) < Q)		_logLik += newLL - oldLL;
			else						_regionExpR[iRegion] = oldValue;
		}
	return 1;
}

//// updates proposal step size, rate for random walk
void updateRate(double &rate, double currAR, double optAR, double delta)
{
	rate = rate * (1.0 + delta * (currAR - optAR));
}

void updateDensityGT() //// updates _densityGT vector (probability density for generation time)
{
	int day;
	double a		= pow(mGT() / sdGT(), 2);
	double b		= mGT() / pow(sdGT(), 2);
	double sumProba = pgamma(_maxDuration, a, b); // pgamma gives the cumulative distribution
	for (day = 0; day < _maxDuration; ++day)
		_densityGT[day] = (pgamma((double)(day + (int)1), a, b) - pgamma(day, a, b)) / sumProba;
}

int drawHospitalInRegion		(int region)
{
	return _hospitalInRegion[region][int(runif()*_numberOfHospitalsPerRegion[region])];
}
int drawHospitalOutsideRegion	(int region)
{
	int hospital = int(runif() * _numberOfHospitals);
	while (_regionOfHospital[hospital] == region) hospital = int(runif()*_numberOfHospitals);
	return hospital;
}
int drawHospitalSecCase			(int hospitalCase, int H2Htype)
{
	if (H2Htype == SameHosptial	)	return hospitalCase;
	if (H2Htype == SameRegion	)	return drawHospitalInRegion(_regionOfHospital[hospitalCase]);
	else							return drawHospitalOutsideRegion(_regionOfHospital[hospitalCase]);
}

void updateProbaOfTimingTrans	() //// i.e. updates _probaOfTimingTrans and _weightedProbaOfTimingTrans for each infectee/infector (i.e. the w_n from page 2 of PNAS supporting info (section "Update of the source of infection of case i"))
{
	//// this function called after various parameters and quantities have been changed (e.g. _clusterExpR changed in updateExpR() function call.)

	double sumWeight;
	int H2Htype;
	for (int iInfectee = 1; iInfectee < _numberOfCases; ++iInfectee)
	{
		sumWeight = 0.;
		//// calculate weights and sum
		for (set<int>::iterator pInfector = _possibleInfectors[iInfectee].begin(); pInfector != _possibleInfectors[iInfectee].end(); ++pInfector)
		{
			if ((*pInfector) == AnimalReservoir)	_probaOfTimingTrans[iInfectee][*pInfector] = _beta;// alpha()/_numberOfHospitals;//_beta;//
			else
			{
				H2Htype = computeH2Htype(iInfectee, *pInfector);
				_probaOfTimingTrans[iInfectee][*pInfector] = expectedRForCase(*pInfector, H2Htype) * coeffHospitalInteraction(iInfectee, (*pInfector)) * densityGT(_onset[iInfectee] - _onset[*pInfector]);
			}
			sumWeight += _probaOfTimingTrans[iInfectee][*pInfector];
		}
		/// normalise weights using sum
		for (set<int>::iterator pInfector = _possibleInfectors[iInfectee].begin(); pInfector != _possibleInfectors[iInfectee].end(); ++pInfector)
			if ((*pInfector) == AnimalReservoir) _weightedProbaOfTimingTrans[iInfectee][*pInfector] = _probaOfTimingTrans[iInfectee][*pInfector] / sumWeight;
			else 
			{
				H2Htype = computeH2Htype(iInfectee, *pInfector);
				_weightedProbaOfTimingTrans[iInfectee][*pInfector] = _probaOfTimingTrans[iInfectee][*pInfector] / sumWeight;
			}
	}
}
void initiateInfector			()
{
	// Assign infector to each infectee
	int iInfector;
	for (int iInfectee = 1; iInfectee < _numberOfCases; ++iInfectee)
	{
		iInfector				= drawInfector(iInfectee);
		_infector[iInfectee]	= iInfector;
		_secondaryCases[iInfector].insert(iInfectee);
		if (iInfector == AnimalReservoir) // Increment number of introductions for that day
			++_numberOfIntroducers[_onset[iInfectee]];
		else				// Increment the appropriate observed reproduction number of that infector
			++_individualObsR[iInfector][computeH2Htype(iInfectee, iInfector)];
	}

	// Assign generation to each infectee
	int generation, currInfectee;
	for (int iInfectee = 1; iInfectee < _numberOfCases; ++iInfectee)
	{
		generation		= 0;
		currInfectee	= iInfectee;
		// We move up the tree until we hit case 0
		while (_infector[currInfectee] != 0)
		{
			++generation;
			currInfectee = _infector[currInfectee];
		}
		_generation[iInfectee] = generation;
	}

	_numberOfIntroducersWith0Sec = 0;
	for (int iInfector = 1; iInfector < _numberOfCases; ++iInfector)
		if (_infector[iInfector] == AnimalReservoir)	// i.e. if infector is an introducer
			if (_individualObsR[iInfector][0] == 0)		// and if the observed reproduction number of that infector is zero.
				_numberOfIntroducersWith0Sec++; 
}
int drawInfector				(int iInfectee)
{
	return rmultinomialOneValue(_weightedProbaOfTimingTrans[iInfectee], _possibleInfectors[iInfectee]);
}

void fixGenerationDownStream	(int iCase, int generation)
{
	//// i.e. assign to iCase a new value of generation, namely the generation of the infector plus one (by definition of generation). 
	_generation[iCase] = generation;
	if (_secondaryCases[iCase].empty()) return; //// if no secondary cases, can stop here...
	//// .... but if not, loop through all secondary cases and amend their generations. And then recursively amend the generations of their secondary cases, and so on until no more secondary cases. 
	for (set<int>::iterator it = _secondaryCases[iCase].begin(); it != _secondaryCases[iCase].end(); ++it)
		fixGenerationDownStream(*it, generation + 1);
}
void changeInfector				(int iCase, int oldInfector, int newInfector)
{
	//// changes the global variables: _infector, _secondaryCases, _individualObsR, _numberOfIntroducers, _generation
	_infector[iCase] = newInfector;
	_secondaryCases[oldInfector].erase(iCase);
	_secondaryCases[newInfector].insert(iCase);
	if (oldInfector > 0)
	{
		int oldH2Htype = computeH2Htype(iCase, oldInfector);
		--_individualObsR[oldInfector][oldH2Htype];
	}
	if (newInfector > 0)
	{
		int newH2Htype = computeH2Htype(iCase, newInfector);
		++_individualObsR[newInfector][newH2Htype];
	}
	if ((oldInfector != AnimalReservoir) && (newInfector == AnimalReservoir))
	{
		++_numberOfIntroducers[_onset[iCase]];
		if (_individualObsR[iCase][0] == 0) _numberOfIntroducersWith0Sec++;
	}
	else if ((oldInfector == AnimalReservoir) && (newInfector != AnimalReservoir))
	{
		--_numberOfIntroducers[_onset[iCase]];
		if (_individualObsR[iCase][0] == 0) _numberOfIntroducersWith0Sec--;
	}
	// Fix generation for all cases down the tree (for tree root _generation[0] = -1)
	if (_trackGeneration == 1) fixGenerationDownStream(iCase, _generation[newInfector] + 1); //// i.e. assign to iCase a new value of generation, namely the generation of the infector plus one (by definition of generation). 
}
double updateInfector			(int iInfectee) // the double that this function returns is less important than the global variables it changes. 
{
	//// structure of this function is as follows: 
	////	i)		draw a potential new infector
	////	ii)		identify and calculate contritubution to likelihood associated with this iInfectee, the current oldInfector and the potential newInfector (assuming old and new are different - no point otherwise). 
	////	iii)	change oldInfector to newInfector using changeInfector, which changes quanties _infector, _secondaryCases, _individualObsR, _numberOfIntroducers, _generation (and all of these quantities downstream).
	////	iv)		identify and calculate contritubution to likelihood associated with this iInfectee, now that newInfector and oldInfector have been swapped. (and all associated quantities have been changed)
	////	v)		Standard ACCEPT/REJECT MCMC step: if ACCEPT increment global _logLik variable (total likelihood over all patients) by new contribution and return 1 acceptance; if REJECT reset quantities using changeInfector function and return 0 acceptances.

	int oldInfector = _infector[iInfectee];
	int newInfector = drawInfector(iInfectee);
	if (oldInfector == newInfector) return 0.; //// i.e. and don't do the rest of this function (but do add zero to numberacceptances etc.)
	
	// Log-likelihood before change
	double oldLogLik = contribLogLikCaseForUpdateInfector(iInfectee, oldInfector, newInfector);
	oldLogLik += ldensityUndetectedIntroducersWith0Sec(); //// body of ldensityUndetectedIntroducersWith0Sec function commented out - just a relic of older analysis. 

	changeInfector(iInfectee, oldInfector, newInfector);

	// Log-likelihood after change
	double newLogLik = contribLogLikCaseForUpdateInfector(iInfectee, newInfector, oldInfector);
	newLogLik += ldensityUndetectedIntroducersWith0Sec();  //// body of ldensityUndetectedIntroducersWith0Sec function commented out - just a relic of older analysis.

	double logPAcceptance	= newLogLik - oldLogLik;
	double logProposal		= log(_weightedProbaOfTimingTrans[iInfectee][oldInfector]) - log(_weightedProbaOfTimingTrans[iInfectee][newInfector]);
	
	if (log(runif()) < (logPAcceptance + logProposal))			///// i.e. if you accept
	{
		_logLik += logPAcceptance;								//// update _logLik
		return 1.;												//// return 1 which			ADDS		to running total of number of acceptances (numberAcceptedUpdateInfector) when this function is called 
	}
	else   														///// i.e. if you reject
	{
		changeInfector(iInfectee, newInfector, oldInfector);	//// Change infector back again
		return 0.;												//// return 0 which		DOES NOT ADD	to running total of number of acceptances (numberAcceptedUpdateInfector) when this function is called 
	}
}

void simulateEpidemic			()
{
	//// clear everything
	for (int day = 0; day<_maxDay; day++) _simulCaseOnDay[day].clear();
	_simulCaseID.		clear();
	_simulOnset.		clear();
	_simulInfectorID.	clear();
	_simulHospID.		clear();
	_simulRegionID.		clear();
	_simulNbSecCases.	clear();

	int currentCaseID = -1;
	int nbDailyIntro;

	// draw all introductions
	for (int day = 0; day < _maxDay; day++)
	{
		if (_introOverdispersed == 0) nbDailyIntro = poidev(alpha(day));
		else 
		{
			double k		= alphaK();
			double p		= 1. / (1. + k / alpha(day));
			nbDailyIntro	= rnegbin(k, p);
		}
		for (int iCase = 0; iCase < nbDailyIntro; iCase++)
		{
			currentCaseID++;
			_simulCaseID.		push_back(currentCaseID);
			_simulOnset.		push_back(day);
			_simulInfectorID.	push_back(0);
			_simulHospID.		push_back(int(runif() * _numberOfHospitals));
			_simulRegionID.		push_back(_regionOfHospital[_simulHospID[currentCaseID]]);
			_simulNbSecCases.	push_back({ 0,0,0 });

			_simulCaseOnDay[day].insert(currentCaseID);
		}
	}

	_simulNbIntro = currentCaseID;
	int nSecCase = 0;
	// draw secondary cases and transmission tree
	for (int day = 0; day < _maxDay; day++)
		if (_simulCaseOnDay[day].size() > 0)
			for (set<int>::iterator iCase = _simulCaseOnDay[day].begin(); iCase != _simulCaseOnDay[day].end(); ++iCase)
				for (int H2Htype = 0; H2Htype<_numberOfH2Htype; H2Htype++)
				{
					if ((H2Htype>0) | (_variabilityRAtClusterLevel == vC_No_Hetero)) _simulNbSecCases[*iCase][H2Htype] = poidev(mR(H2Htype));
					else 
					{
						double p = 1. / (1. + kR(H2Htype) / mR(H2Htype));
						_simulNbSecCases[*iCase][H2Htype] = rnegbin(kR(H2Htype), p);
					}
					nSecCase += _simulNbSecCases[*iCase][H2Htype];

					for (int iSecCase = 0; iSecCase<_simulNbSecCases[*iCase][H2Htype]; iSecCase++)
					{
						currentCaseID++;
						_simulCaseID.		push_back(currentCaseID);

						int dayOnsetSecCase = day + rmultinomialOneValue(_densityGT, _maxDuration);
						
						_simulOnset.		push_back(dayOnsetSecCase);
						_simulInfectorID.	push_back(*iCase);
						_simulHospID.		push_back(drawHospitalSecCase(_simulHospID[*iCase], H2Htype));
						_simulRegionID.		push_back(_regionOfHospital[_simulHospID[currentCaseID]]);
						_simulNbSecCases.	push_back({ 0,0,0 });
						if (dayOnsetSecCase<_maxDay) _simulCaseOnDay[dayOnsetSecCase].insert(currentCaseID);
					}
				}
	_simulNbCase = currentCaseID + 1;
}
void simulateEpidemicWithCluster()
{
	int maxOutbreakSize			= 50000;
	int maxOutbreakSizeReached	= 0;

	//// clear everything
	for (int day = 0; day < _maxDay; day++) _simulCaseOnDay[day].	clear();
	_simulCaseID.				clear();
	_simulOnset.				clear();
	_simulInfectorID.			clear();
	_simulHospID.				clear();
	_simulRegionID.				clear();
	_simulNbSecCases.			clear();
	_simulClusterOfCase.		clear();
	_simulClusterR.				clear();
	_simulClusterStart.			clear();
	_simulClusterHospital.		clear();
	_simulLastOnsetInHospital.	clear();
	_simulLastClusterInHospital.clear();
	_simulNbInClusterBeforeDay.	clear();
	_simulNbInClusterOnDay.		clear();

	int currentCaseID = -1;
	int nbDailyIntro;

	// draw all introductions
	for (int day = 0; day < _maxDay; day++)
	{
		if (_introOverdispersed == 0) nbDailyIntro = poidev(alpha(day));
		else 
		{
			double k = alphaK();
			double p = 1. / (1. + k / alpha(day));
			nbDailyIntro = rnegbin(k, p);
		}
		for (int iCase = 0; iCase < nbDailyIntro; iCase++)
		{
			currentCaseID++;
			_simulCaseID.push_back(currentCaseID);
			_simulOnset.push_back(day);
			_simulInfectorID.push_back(0);
			_simulHospID.push_back(int(runif()*_numberOfHospitals));
			_simulRegionID.push_back(_regionOfHospital[_simulHospID[currentCaseID]]);
			_simulNbSecCases.push_back({ 0,0,0 });
			_simulClusterOfCase.push_back(-1);

			_simulCaseOnDay[day].insert(currentCaseID);
		}
	}

	for (int iRegion = 0; iRegion < _numberOfRegions; iRegion++)
		if (_heterogeneityRegionR == 1)	_simulRegionR[iRegion] = _regionExpR[iRegion];		else	_simulRegionR[iRegion] = mR(1);

	for (int hospitalID = 0; hospitalID < _numberOfHospitals; hospitalID++)
	{
		_simulLastOnsetInHospital.	push_back(-1000);
		_simulLastClusterInHospital.push_back(-1);
	}

	_simulNbIntro = currentCaseID;
	// draw secondary cases and transmission tree
	int H2Htype;
	int currentClusterID = -1;
	for (int day = 0; day < _maxDay; day++)
	{
		for (int iCluster = 0; iCluster<currentClusterID; iCluster++)
		{
			int delaySinceLastCase = day - _simulLastOnsetInHospital[_simulClusterHospital[iCluster]];
			if (delaySinceLastCase<_minDelayBetweenClusters) _simulNbInClusterOnDay[iCluster] = 0;
		}

		if (_simulCaseOnDay[day].size()>0)
		{
			for (set<int>::iterator iCase = _simulCaseOnDay[day].begin(); iCase != _simulCaseOnDay[day].end(); ++iCase)
			{
				int delaySinceLastCase = day - _simulLastOnsetInHospital[_simulHospID[*iCase]];
				if (delaySinceLastCase >= _minDelayBetweenClusters)
				{
					currentClusterID++;

					_simulClusterStart.push_back(day);
					_simulClusterHospital.push_back(_simulHospID[*iCase]);

					_simulLastOnsetInHospital[_simulHospID[*iCase]] = day;
					_simulLastClusterInHospital[_simulHospID[*iCase]] = currentClusterID;

					_simulClusterR.push_back({ 0,0,0 });
					H2Htype = SameHosptial;
					if (_heterogeneityR == 0) _simulClusterR[currentClusterID][H2Htype] = mR(H2Htype);
					else {
						double p = 1. / (1. + kR(H2Htype) / mR(H2Htype));
						_simulClusterR[currentClusterID][H2Htype] = rgamma(aR(H2Htype), bR(H2Htype));
					}
					_simulClusterR[currentClusterID][1] = mR(1);
					_simulClusterR[currentClusterID][2] = mR(2);
					_simulNbInClusterBeforeDay.push_back(0);
					_simulNbInClusterOnDay.push_back(0);
				}

				_simulClusterOfCase[*iCase] = _simulLastClusterInHospital[_simulHospID[*iCase]];
				_simulNbInClusterOnDay[_simulClusterOfCase[*iCase]]++;

				H2Htype = SameHosptial;
				_simulNbSecCases[*iCase][H2Htype]	= poidev(_simulClusterR[_simulClusterOfCase[*iCase]][H2Htype] * reductionRForCaseWithCumNbCases(_simulNbInClusterBeforeDay[_simulClusterOfCase[*iCase]]));
				_simulNbSecCases[*iCase][1]			= poidev(_simulRegionR[_simulRegionID[*iCase]]);
				_simulNbSecCases[*iCase][2]			= poidev(mR(2));

				for (H2Htype = 0; H2Htype < _numberOfH2Htype; H2Htype++)
				{
					for (int iSecCase = 0; iSecCase < _simulNbSecCases[*iCase][H2Htype]; iSecCase++)
					{
						currentCaseID++;
						if (currentCaseID>maxOutbreakSize)
						{
							maxOutbreakSizeReached = 1;
							break;
						}
						if (currentCaseID)
							_simulCaseID.push_back(currentCaseID);

						int dayOnsetSecCase = day + rmultinomialOneValue(_densityGT, _maxDuration);
						_simulOnset.push_back(dayOnsetSecCase);

						_simulInfectorID.push_back(*iCase);
						_simulHospID.push_back(drawHospitalSecCase(_simulHospID[*iCase], H2Htype));
						_simulRegionID.push_back(_regionOfHospital[_simulHospID[currentCaseID]]);

						_simulClusterOfCase.push_back(-1);

						_simulNbSecCases.push_back({ 0,0,0 });
						if (dayOnsetSecCase<_maxDay) _simulCaseOnDay[dayOnsetSecCase].insert(currentCaseID);
					}
					if (maxOutbreakSizeReached == 1) break;
				}
				if (maxOutbreakSizeReached == 1) break;
			}

			for (int iCluster = 0; iCluster < currentClusterID; iCluster++)
			{
				int delaySinceLastCase = day - _simulLastOnsetInHospital[_simulClusterHospital[iCluster]];
				if (delaySinceLastCase<_minDelayBetweenClusters) _simulNbInClusterBeforeDay[iCluster] += _simulNbInClusterOnDay[iCluster];
			}
		}
		if (maxOutbreakSizeReached == 1) break;
	}

	_simulNbCase			= currentCaseID + 1;
	_simulNbCluster			= currentClusterID + 1;

	_simulMeanClusterSize	= 0;
	_simulMaxClusterSize	= 0;
	_simulPClusterSize1		= 0;
	_simulPClusterSize10	= 0;
	for (int iCluster = 0; iCluster<_simulNbCluster; iCluster++)
	{
		_simulMeanClusterSize	+= _simulNbInClusterBeforeDay[iCluster];
		_simulPClusterSize1		+= (_simulNbInClusterBeforeDay[iCluster] == 1);
		_simulPClusterSize10	+= (_simulNbInClusterBeforeDay[iCluster] > 10);
		if (_simulNbInClusterBeforeDay[iCluster]>_simulMaxClusterSize) _simulMaxClusterSize = _simulNbInClusterBeforeDay[iCluster];
	}
	_simulMeanClusterSize	= _simulMeanClusterSize / _simulNbCluster;
	_simulPClusterSize1		= _simulPClusterSize1	/ _simulNbCluster;
	_simulPClusterSize10	= _simulPClusterSize10	/ _simulNbCluster;

	int durationPeriod				= 60;
	_simulMaxNbCaseOverPeriod		= 0;
	_simulMaxNbClusterOverPeriod	= 0;

	int iSubsequentCluster;
	int nbCaseInPeriod, nbClusterInPeriod;
	for (int iCluster = 0; iCluster < _simulNbCluster; iCluster++)
	{
		nbCaseInPeriod		= 0;
		nbClusterInPeriod	= 0;
		for (iSubsequentCluster = iCluster; iSubsequentCluster < _simulNbCluster; iSubsequentCluster++)
			if (_simulClusterStart[iSubsequentCluster] <= _simulClusterStart[iCluster] + durationPeriod)
			{
				nbCaseInPeriod += _simulNbInClusterBeforeDay[iSubsequentCluster];
				nbClusterInPeriod++;
			}
			else 
			{
				if (_simulMaxNbCaseOverPeriod<nbCaseInPeriod) _simulMaxNbCaseOverPeriod = nbCaseInPeriod;
				if (_simulMaxNbClusterOverPeriod<nbClusterInPeriod) _simulMaxNbClusterOverPeriod = nbClusterInPeriod;
				break;
			}
	}
}
