#include "Libraries.h"
#include "Macros.h"
#include "DataStructure.h"
#include "randlib_par.h"
#include "DJL_Structs.h"
#include "Globals.h"
#include "Likelihood.h"

// =============================================================================
// Common functions implementation
// =============================================================================

double mR(int H2Htype) //// returns the current values of reproduction numbers R_c, R_R or R_O, depending on H2Htype
{
	return _parameter[(int)2 + H2Htype];
}
double kR(int H2Htype) // Dispersion parameter
{
	if (_cvParam == 0)	return _parameter[(int)5 + H2Htype];
	else				return 1. / pow(_parameter[(int)5 + H2Htype], 2);
}
double aR(int H2Htype)
{
	return kR(H2Htype);	//option k
}
double bR(int H2Htype)
{
	return kR(H2Htype) / mR(H2Htype);	//option k
}
double mReductionR()
{
	return _parameter[10];
}
double kReductionR()
{
	return 1./pow(_parameter[11],2);
}
double aReductionR()
{
	return kReductionR();
}
double bReductionR()
{
	return kReductionR()/mReductionR();
}
double reductionRForCaseWithCumNbCases(int nbCumNbCases)
{
	double baselineNb = 10;
	double reductionAfterBaselineNb;
	reductionAfterBaselineNb = mReductionR();
	return pow(1 + nbCumNbCases, -(-log(reductionAfterBaselineNb) / log(1 + baselineNb)));
}
double reductionRForDelay(int delay)
{
	double baselineDelay = 14;
	double reductionAfterBaselineDelay;
	reductionAfterBaselineDelay = mReductionR();
	return pow(1 + delay, -(-log(reductionAfterBaselineDelay) / log(1 + baselineDelay)));
}
double reductionRForCase(int iCase)
{
		 if (_withReductionR == 0)	return 1.;
	else if (_withReductionR == 1)	return reductionRForCaseWithCumNbCases(_rankOfCaseInCluster[iCase]);
	else							return reductionRForDelay(_onset[iCase] - _startOfCluster[_clusterOfCase[iCase]]);
}
double expectedRForCase(int iCase, int H2Htype)
{
	if (H2Htype == SameHosptial)
	{
			 if (_variabilityRAtClusterLevel == vC_No_Hetero		)	return mR(H2Htype) * reductionRForCase(iCase);
		else if (_variabilityRAtClusterLevel == vC_IndividualLevel	)	return _individualExpR[iCase][H2Htype]				* reductionRForCase(iCase);
		else															return _clusterExpR[_clusterOfCase[iCase]][H2Htype] * reductionRForCase(iCase);
	}
	else 
	{
		if ((H2Htype == SameRegion) & (_heterogeneityRegionR == 1)) return _regionExpR[_region[iCase]];
		else return mR(H2Htype);
	}
}
double alpha(int day)
{
	double effect = _parameter[7] * exp(_parameter[8] * day);
	if (_withSeasonality == No_Seasonality) return effect;
	else if (_withSeasonality == Seasonality_AboveOne)
	{
		double time				= double(day) / 365;
		double effectSeasonal	= (1 + _parameter[13] * cos(2 * 3.141592654*(time - _parameter[12])));
		if (effectSeasonal < 1) effectSeasonal = 1;
		return effect * effectSeasonal;
	}
	else if (_withSeasonality == Seasonality_TrueCos)
	{
		double time				= double(day) / 365;
		double effectSeasonal	= (1 + _parameter[13] * cos(2 * 3.141592654*(time - _parameter[12])));
		effect					= effect * effectSeasonal;
		if (effect < 0.000000001) effect = 0.000000001;
		return effect;
	}
	else 
	{
		double time				= double(day) / 365;
		double effectSeasonal	= _parameter[13] * cos(2 * 3.141592654*(time - _parameter[12]));
		if (effectSeasonal < 0) effectSeasonal = 0;
		effect					= effect + effectSeasonal;
		return effect;
	}
	return effect;
}
double alphaK()
{
	return _parameter[9];
}
double mGT() // Mean generation time
{
	return _parameter[0];
}
double sdGT() // Generation time std
{
	return _parameter[1];
}
double densityGT(int x)
{
	if (x < 0)				return 0.;
	if (x >= _maxDuration)	return 0.;
	return _densityGT[x];
}


void initiateVariabilityExpR() //// changes _individualExpR OR _clusterExpR (SameHospital) (at cluster level, as opposed to region or hospital) for each individual / cluster depending on "scenario variable" _variabilityRAtClusterLevel  
{
	int H2Htype = SameHosptial;
	double sumObsRCluster, nbCasesCluster;
	if (_variabilityRAtClusterLevel == vC_IndividualLevel) // individual level
		for (int iInfectee = 1; iInfectee < _numberOfCases; ++iInfectee)
			_individualExpR[iInfectee][H2Htype] = rgamma(aR(H2Htype) + _individualObsR[iInfectee][H2Htype], bR(H2Htype) + 1);
	else if (_variabilityRAtClusterLevel == vC_ClusterLevel) // cluster level
	{
		if (_heterogeneityR == 1)
			for (int iCluster = 1; iCluster < _numberOfClusters; ++iCluster)
			{
				sumObsRCluster = 0;
				nbCasesCluster = 0;
				for (set<int>::iterator it = _casesInCluster[iCluster].begin(); it != _casesInCluster[iCluster].end(); ++it)
				{
					nbCasesCluster++;
					sumObsRCluster += _individualObsR[*it][H2Htype];
				}
				_clusterExpR[iCluster][H2Htype] = rgamma(aR(H2Htype) + sumObsRCluster, bR(H2Htype) + nbCasesCluster);
			}
		if (_heterogeneityReductionR == 1)
			for (int iCluster = 1; iCluster < _numberOfClusters; ++iCluster)
				_clusterReductionR[iCluster] = rgamma(aReductionR(), bReductionR());
	}
	H2Htype = SameRegion;
	if (_heterogeneityRegionR == 1)
		for (int iRegion = 0; iRegion < _numberOfRegions; iRegion++)
			_regionExpR[iRegion] = rgamma(aR(H2Htype), bR(H2Htype));
}

/////////////////		/////////////////		/////////////////		/////////////////		/////////////////		
/////////////////		Log likelihood functions						/////////////////

double probaOfTimingTrans					(int iInfectee, int iInfector)
{
	return coeffHospitalInteraction(iInfectee, iInfector) * densityGT(_onset[iInfectee] - _onset[iInfector]);
}
double ldensityUndetectedIntroducersWith0Sec()
{
	return 0;
}
double ldensityIndividualR					(int iCase, int H2Htype)
{
	return  dlpois(_individualObsR[iCase][H2Htype], expectedRForCase(iCase, H2Htype));
}
double logProbaIntro						(int day)
{
	double LL;
	if (_introOverdispersed == 0) LL = dlpois(_numberOfIntroducers[day], alpha(day));
	else
	{
		double k = alphaK();
		double p = 1. / (1. + k / alpha(day));
		LL = dlnegbin(_numberOfIntroducers[day], k, p);
	}
	if (_numberOfIntroducers[day] > 0)
	{
		if (_spatialLevel >= SL_Hospital) LL += -_numberOfIntroducers[day] * log(_numberOfHospitals);
		else if (_spatialLevel == SL_Region) LL += -_numberOfIntroducers[day] * log(_numberOfRegions);
	}
	return LL;
}
double logLik()
{
	int H2Htype		= SameHosptial;
	double logLik	= 0;

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** 
	//// **** L^cluster (from page 1 PNAS supporting info)
	if (_variabilityRAtClusterLevel == vC_IndividualLevel)		// individual
		for (int iInfector = 1; iInfector < _numberOfCases; ++iInfector)
			logLik += ldgamma(_individualExpR[iInfector][H2Htype], aR(H2Htype), bR(H2Htype));
	else if (_variabilityRAtClusterLevel == vC_ClusterLevel)	// cluster
	{
		if (_heterogeneityR == 1)
			for (int iCluster = 0; iCluster < _numberOfClusters; ++iCluster)
				logLik += ldgamma(_clusterExpR[iCluster][H2Htype]	, aR(H2Htype)	, bR(H2Htype));
		if (_heterogeneityReductionR == 1)
			for (int iCluster = 0; iCluster < _numberOfClusters; ++iCluster)
				logLik += ldgamma(_clusterReductionR[iCluster]		, aReductionR()	, bReductionR());
	}

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** 
	//// **** L_n^trans  (from page 1 supporting info)
	H2Htype = SameRegion;
	if (_heterogeneityRegionR == 1)
		for (int iRegion = 0; iRegion<_numberOfRegions; iRegion++)
			logLik += ldgamma(_regionExpR[iRegion], aR(H2Htype), bR(H2Htype));

	//// first line of L_n^trans definition (summed / multiplied over case numbers n and H2Htypes). 
	for (int iInfector = 1; iInfector < _numberOfCases; ++iInfector) //logLik += lProbaDetectionIntroduction(iInfector);
		for (H2Htype = 0; H2Htype<_numberOfH2Htype; H2Htype++)
			logLik += ldensityIndividualR(iInfector, H2Htype);

	//// second and third lines of L_n^trans definition (summed / multiplied over case numbers n). Serial interval of timing of secondary infections from infectee, and uniform distribution of hospital interactions. 
	for (int iInfectee = 1; iInfectee < _numberOfCases; ++iInfectee)
		if (_infector[iInfectee] > 0) logLik += log(probaOfTimingTrans(iInfectee, _infector[iInfectee]));
	
	//// body of code of this function commented out - just returns zero and no mention made of it in paper - relic of older analysis. 
	logLik += ldensityUndetectedIntroducersWith0Sec(); 

	//// **** //// **** //// **** //// **** //// **** //// **** //// **** 
	//// **** prod_t( L_t^intro )  (from page 1 supporting info)
	for (int day = _minDay; day <= _maxDay; ++day)	logLik += logProbaIntro(day);

	return logLik;
}
void contribLogLikDownStream				(int iCase, double &contrib)
{
	//// Function calculates contribution of case n to L_n^trans part of likelihood (summed / multiplied over H2Htypes). Effectively returns a double. 

	for (int H2Htype = 0; H2Htype < _numberOfH2Htype; H2Htype++) contrib += ldensityIndividualR(iCase, H2Htype); //// first line of L_n^trans definition (summed / multiplied over H2Htypes). 
	if (_trackGeneration == 1)
	{
		if (_secondaryCases[iCase].empty()) return;  //// if there are no secondary cases, can stop here...
		//// .... but if not, loop through all secondary cases and amend their contributions to the likelihood. And then recursively amend the generations of their secondary cases, and so on until no more secondary cases. 
		for (set<int>::iterator it = _secondaryCases[iCase].begin(); it != _secondaryCases[iCase].end(); ++it)
			contribLogLikDownStream(*it, contrib);
	}
}
double contribLogLikCaseForUpdateInfector	(int iCase, int oldInfector, int newInfector)
{
	// Contribution of case iCase to log-likelihood, given oldInfector and newInfector. Not to be confused with total log likelihood variable _logLik (which is global)
	double logLik_iCase = 0.;
	contribLogLikDownStream(iCase, logLik_iCase); //// Function calculates contribution of case n to L_n^trans part of likelihood (summed / multiplied over H2Htypes), as well as all downstream cases affected.
	if (oldInfector > 0)
	{
		for (int H2Htype = 0; H2Htype < _numberOfH2Htype; H2Htype++)
			logLik_iCase += ldensityIndividualR(oldInfector, H2Htype);
		logLik_iCase += log(probaOfTimingTrans(iCase, oldInfector)); //// (log) second and third lines of L_n^trans definition. Serial interval of timing of secondary infections from infectee, and uniform distribution of hospital interactions. 
	}
	if (newInfector > 0)
	{	
		for (int H2Htype = 0; H2Htype < _numberOfH2Htype; H2Htype++)
			logLik_iCase += ldensityIndividualR(newInfector, H2Htype);
	}
	logLik_iCase += logProbaIntro(_onset[iCase]); //// **** L_t^intro ( from page 1 supporting info)
	return logLik_iCase;
}

int computeH2Htype(int case1, int case2) // CAREFUL: CHANGE CoeffHospInter and logProbaIntro if this function changes
{
	//_spatialLevel: 0: no space / 1: region / 2: hospital / 3: region+hospital
	if (_spatialLevel == SL_No_Space)  // 0: no space
		return SameHosptial;

	if (_spatialLevel == SL_Region)
	{
		if (_region[case1] == _region[case2]) return SameHosptial;
		else return SameRegion;
	}

	if (_spatialLevel == SL_Hospital)
	{
		if (_hospital[case1]	== _hospital[case2])	return SameHosptial;
		else return SameRegion;
	}

	if (_spatialLevel == SL_RegionAndHospital)
	{
		if (_hospital[case1]	== _hospital[case2])	return SameHosptial;
		if (_region[case1]		== _region[case2])		return SameRegion	;
		return DifferentRegion; // other region
	}
	return -1;
}
double coeffHospitalInteraction(int iInfectee, int iInfector)
{
	int H2Htype = computeH2Htype(iInfectee, iInfector);

	if (H2Htype == SameHosptial) return 1;

	if (H2Htype == SameRegion)
	{
		if (_spatialLevel == SL_Region	) return 1. / double(_numberOfRegions		- 1.);
		if (_spatialLevel == SL_Hospital) return 1. / double(_numberOfHospitals		- 1.);
		return	1. / double(_numberOfHospitalsPerRegion[_region[iInfectee]]	- 1.);
	}

	if (H2Htype == DifferentRegion)			return 1. / (double)((double)_numberOfHospitals - (double)_numberOfHospitalsPerRegion[_region[iInfectee]]); //spatial level 3
	return -1;
}

