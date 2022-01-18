#pragma once

#include "Libraries.h"

double updateParameterNormalProposal(int parameterNumber, double (*logLik)());
double updateParameterLogNormalProposal(int parameterNumber, double (*logLik)()); 
double updateParameterGT(int parameterNumber, double(*logLik)()); 

//// Depending on the parameter number, this calls either of the three functions above. 
double updateParameter(int parameterNumber); 

//// updateExpR not an paramter update function per se, as it assumes parameters have already gone through an ACCEPT/REJECT step.
///// Instead, updateExpR updates varrious quantities associated with the likelihood and that are necessary to calculate the likelihood, e.g. the expected R's at cluster level (or individual level depending on _variabilityRAtClusterLevel)
double updateExpR();

//// updates proposal step size, rate for random walk
void updateRate(double& rate, double currAR, double optAR, double delta);
void updateDensityGT(); //// updates _densityGT vector (probability density for generation time)

int drawHospitalSecCase(int hospitalCase, int H2Htype);

void	updateProbaOfTimingTrans();
void initiateInfector			();
int	drawInfector(int iInfectee);

void fixGenerationDownStream	(int iCase, int generation);
void changeInfector				(int iCase, int oldInfector, int newInfector);
double updateInfector			(int iInfectee) ; // the double this function returns is less important than the global variables it changes. 

void simulateEpidemic(); 
void simulateEpidemicWithCluster(); 