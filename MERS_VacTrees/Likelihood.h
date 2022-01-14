#pragma once

// =============================================================================
// Common functions declaration
// =============================================================================

double mR(int H2Htype); //// returns the current values of reproduction numbers R_c, R_R or R_O, depending on H2Htype
double kR(int H2Htype);  // Dispersion parameter
double aR(int H2Htype); 
double bR(int H2Htype); 
double mReductionR()  ; 
double kReductionR()  ; 
double aReductionR()  ; 
double bReductionR()  ; 

double reductionRForCaseWithCumNbCases(int nbCumNbCases); 
double reductionRForDelay(int delay); 
double reductionRForCase(int iCase);
double expectedRForCase(int iCase, int H2Htype);

double alpha(int day); 
double alphaK();
double mGT();// Mean generation time
double sdGT(); // Generation time std
double densityGT(int x);

double	mGT(); // Mean generation time
double	sdGT(); // Generation time std
void	updateDensityGT();
double	densityGT(int x);

int computeH2Htype(int case1, int case2); // CAREFUL: CHANGE CoeffHospInter and logProbaIntro if this function changes
double coeffHospitalInteraction(int iInfectee, int iInfector);


double probaOfTimingTrans(int iInfectee, int iInfector); 
double ldensityUndetectedIntroducersWith0Sec()										  ; 
double ldensityIndividualR(int iCase, int H2Htype)									  ; 
double logProbaIntro(int day)														  ; 
double logLik()																		  ; 
void contribLogLikDownStream(int iCase, double& contrib)							  ; 
double contribLogLikCaseForUpdateInfector(int iCase, int oldInfector, int newInfector); 

// Move this to Initialize functions.
void initiateVariabilityExpR(); 

