#pragma once

#include "Libraries.h"
#include "Structs.h"


bool FileExists(string FileName);
void ReadInParams(ModelRun& MR, string pParamFileName);

VaccCampaignStrategy	Convert_VaccCampaignStrategy_FromString	(const std::string& OptionString); 
ReactiveLevel			Convert_ReactiveLevel_FromString		(const std::string& OptionString);
std::string Convert_VaccCampaignStrategy_FromEnumClass	(VaccCampaignStrategy VacCampStrategy);
std::string Convert_ReactiveLevel_FromEnumClass			(ReactiveLevel ReactLevel);

void PutColNamesOn_Chains		(ofstream& Chains);
void PutColNamesOn_IndividualR(ofstream& IndividualR);
void PutColNamesOn_Trees		(ofstream& Trees);
void PutColNamesOn_Trees_CF		(ofstream& outputFileTree_CF);
void PutColNamesOn_SimulCluster	(ofstream& outputFileSimulCluster);
void PutColNamesOn_AllOutput	(AllOutput& OUTPUT);
void InitializeOutputValues		(AllOutput& OUTPUT);

void loadAndInitializeData(string fileName, ModelRun& MR, bool printData);

void WriteModelMetaData(AllOutput& OUTPUT, ModelRun& MR, FileStrings_Struct FileStrings, double delta, double kRinit, double pDetectIndex, double kIntro);  //// unseen "arguments" are various global variables, e.g. _spatialLevel

template <typename T>
std::string ToStringWithPrecision(const T a_value, const int n); 
std::string ChooseScenarioName(ModelRun& MR);
int Choose_numberOfHospitals(int dataSevereCases); 
