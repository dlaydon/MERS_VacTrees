#pragma once

using namespace std;

#include "Macros.h"

//// **** //// **** //// **** //// **** //// **** //// **** 
//// **** //// **** //// **** TEMPLATES FOR MEMORY ALLOCATION/POPULATION/CHECKING 


template <typename TYPE> void Allocate_2D_Array		(TYPE **    &OBJECT, int Dim1, int Dim2)
{
	OBJECT = new TYPE *[Dim1]();
	for (int row = 0; row < Dim1; row++)	OBJECT[row] = new TYPE[Dim2]();
}
template <typename TYPE> void Allocate_3D_Array		(TYPE ***   &OBJECT, int Dim1, int Dim2, int Dim3)
{
	OBJECT = new TYPE **[Dim1]();
	for (int row = 0; row < Dim1; row++)	Allocate_2D_Array(OBJECT[row], Dim2, Dim3);
}
template <typename TYPE> void Allocate_4D_Array		(TYPE ****  &OBJECT, int Dim1, int Dim2, int Dim3, int Dim4)
{
	OBJECT = new TYPE ***[Dim1]();
	for (int row = 0; row < Dim1; row++)	Allocate_3D_Array(OBJECT[row], Dim2, Dim3, Dim4);
}
template <typename TYPE> void Allocate_5D_Array		(TYPE ***** &OBJECT, int Dim1, int Dim2, int Dim3, int Dim4, int Dim5)
{
	OBJECT = new TYPE ****[Dim1]();
	for (int row = 0; row < Dim1; row++)	Allocate_4D_Array(OBJECT[row], Dim2, Dim3, Dim4, Dim5);
}

template <typename T> bool IsEqual_1D_Arrays	(T*    const &lhs, T*    const &rhs, int dim)
{
	bool All_Equal = 1;
	if (lhs != 0 && rhs != 0)	All_Equal = (memcmp(lhs, rhs, sizeof(T) * dim) == 0);		//// if both arrays non-empty compare them
	else if (lhs == 0 && rhs == 0)	All_Equal = 1;											//// if both arrays empty return true
	else							All_Equal = 0; 											//// if one array empty and other not then return false. 
	return All_Equal;
}
template <typename T> bool IsEqual_2D_Arrays	(T**   const &lhs, T**   const &rhs, int dim1, int dim2)
{
	bool All_Equal = 1; 
	if (lhs != 0 && rhs != 0)														//// if both arrays non-empty compare them
	{
		for (int row = 0; row < dim1; row++)
			if (IsEqual_1D_Arrays(lhs[row], rhs[row], dim2) == false) { All_Equal = 0; break; }
	}
	else if (lhs == 0 && rhs == 0)
	{
		All_Equal = 1; 																//// if both arrays empty return true
	}
	else All_Equal = 0;																//// if one array empty and other not then return false. 

	return All_Equal; 
}
template <typename T> bool IsEqual_3D_Arrays	(T***  const &lhs, T***  const &rhs, int dim1, int dim2, int dim3)
{
	bool All_Equal = 1;
	if (lhs != 0 && rhs != 0)
	{
		for (int row = 0; row < dim1; row++)
			if (IsEqual_2D_Arrays(lhs[row], rhs[row], dim2, dim3) == false) { All_Equal = 0; break; }
	}
	else if (lhs == 0 && rhs == 0)
	{
		All_Equal = 1; 																//// if both arrays empty return true
	}
	else All_Equal = 0;																//// if one array empty and other not then return false.
	return All_Equal;
}
template <typename T> bool IsEqual_4D_Arrays	(T**** const &lhs, T**** const &rhs, int dim1, int dim2, int dim3, int dim4)
{
	bool All_Equal = 1;
	if (lhs != 0 && rhs != 0)
	{
		for (int row = 0; row < dim1; row++)
			if (IsEqual_3D_Arrays(lhs[row], rhs[row], dim2, dim3, dim4) == false) { All_Equal = 0; break; }
	}
	else if (lhs == 0 && rhs == 0)
	{
		All_Equal = 1; 																//// if both arrays empty return true
	}
	else All_Equal = 0;																//// if one array empty and other not then return false.
	return All_Equal;
}

template <typename T> void SetEqual_1D_Arrays	(T*    &ArrayToChange, T*    const &ArrayToCopy, int dim) //// assumes ArrayToChange and ArrayToCopy have same dimensions. 
{
	if (ArrayToChange != 0) if (ArrayToCopy != 0) for (int element = 0; element < dim; element++) ArrayToChange[element] = ArrayToCopy[element];
}
template <typename T> void SetEqual_2D_Arrays	(T**   &ArrayToChange, T**   const &ArrayToCopy, int dim1, int dim2) //// assumes ArrayToChange and ArrayToCopy have same dimensions. 
{
	if (ArrayToChange != 0) if (ArrayToCopy != 0) for (int element = 0; element < dim1; element++) SetEqual_1D_Arrays(ArrayToChange[element], ArrayToCopy[element], dim2);
}
template <typename T> void SetEqual_3D_Arrays	(T***  &ArrayToChange, T***  const &ArrayToCopy, int dim1, int dim2, int dim3)  //// assumes ArrayToChange and ArrayToCopy have same dimensions. 
{
	if (ArrayToChange != 0) if (ArrayToCopy != 0) for (int element = 0; element < dim1; element++) SetEqual_2D_Arrays(ArrayToChange[element], ArrayToCopy[element], dim2, dim3);
}
template <typename T> void SetEqual_4D_Arrays	(T**** &ArrayToChange, T**** const &ArrayToCopy, int dim1, int dim2, int dim3, int dim4)  //// assumes ArrayToChange and ArrayToCopy have same dimensions. 
{
	if (ArrayToChange != 0) if (ArrayToCopy != 0) for (int element = 0; element < dim1; element++) SetEqual_3D_Arrays(ArrayToChange[element], ArrayToCopy[element], dim2, dim3, dim4);
}

template <typename TYPE> void Populate_1D_Array(TYPE *		&OBJECT, TYPE Value, int Dim1)
{
	for (int row = 0; row < Dim1; row++)	OBJECT[row] = Value;
}
template <typename TYPE> void Populate_2D_Array(TYPE **		&OBJECT, TYPE Value, int Dim1, int Dim2)
{
	for (int row = 0; row < Dim1; row++) Populate_1D_Array(OBJECT[row], Value, Dim2);
}
template <typename TYPE> void Populate_3D_Array(TYPE ***	&OBJECT, TYPE Value, int Dim1, int Dim2, int Dim3)
{
	for (int row = 0; row < Dim1; row++) Populate_2D_Array(OBJECT[row], Value, Dim2, Dim3);
}



struct QuantIndices_Struct {

	int HugeQuantIndex = 40000; //// Do not set to NULL as it evaluates to zero for int's.

	/// _QI suffix refers to "Quantity Index"
	int PropAverted				= HugeQuantIndex;		//// Proportion of cases averted index
	int PropAverted_HCW			= HugeQuantIndex;		//// Proportion of cases in healthcare workers averted index
	int PropAverted_nHCW		= HugeQuantIndex;		//// Proportion of cases in non-healthcare workers averted index
	int DeathsAverted_CF		= HugeQuantIndex;		////  Number of deaths averted index
	int DeathsAverted_CF_HCW	= HugeQuantIndex;		////  Number of deaths in healthcare workers averted index
	int DeathsAverted_CF_nHCW	= HugeQuantIndex;		////  Number of deaths in non-healthcare workers averted index

	int FinalCaseDate_CF		= HugeQuantIndex;		
};
struct SingleOutputStruct {

	string FilePath, FileName, delimiter = "\t";
	int NumRows, NumCols; 
	ofstream Stream;

	bool HasRowNames = false, HasColNames = true;
	std::vector<string> RowNames, Colnames;

	SingleOutputStruct(string FilePathArg, string FileNameArg, int NRows, int NCols)
	{
		FilePath = FilePathArg; 
		FileName = FileNameArg; 
		NumRows = NRows; 
		NumCols = NCols; 
		ofstream Stream(FilePathArg.c_str());
	}
	template <typename Type> void Append(Type ValueToAddOn)
	{
		Stream << ValueToAddOn << delimiter;
	}
	void EndLine()
	{
		Stream << std::endl;
	}

	//// Write output functions. Helpful if all quantities contained in one object, e.g. a matrix of parameter chains. Won't work if you need to write multiple objects to the same stream, e.g. Chains adding parameter values, but also _probaSource and values from  simulation. 
	void WriteOutput(std::vector<vector<double>> TwoDim_Vec) //// version for vector of vectors
	{
		if (HasColNames) for (int col = 0; col < NumCols; col++) Append(Colnames[col]);
		for (int row = 0; row < NumRows; row++)
		{
			if (HasRowNames) Append(RowNames[row]); 
			for (int col = 0; col < NumCols; col++) Append(TwoDim_Vec[row][col]); 
			EndLine(); 
		}
	}
	void WriteOutput(double **TwoDim_Vec)	 //// version for array of arrays
	{
		if (HasColNames) for (int col = 0; col < NumCols; col++) Append(Colnames[col]);
		for (int row = 0; row < NumRows; row++)
		{
			if (HasRowNames) Append(RowNames[row]); 
			for (int col = 0; col < NumCols; col++) Append(TwoDim_Vec[row][col]);
			EndLine(); 
		}
	}
	template <typename Type> void WriteOutput(Type TwoDimContainer)//// template version
	{
		if (HasColNames) for (int col = 0; col < NumCols; col++) Append(Colnames[col]);
		for (int row = 0; row < NumRows; row++)
		{
			if (HasRowNames) Append(RowNames[row]);
			for (int col = 0; col < NumCols; col++) Append(TwoDimContainer[row][col]);
			EndLine();
		}
	}
};
struct FileStrings_Struct {

	string pathInput, pathOutput;
	string inputFileName, Simon_scenarioName, DJL_scenarioName;
	string inputFile;
	string Chains, IndividualR, Infector, ClusterR, ClusterReductionR, SimulCluster, Tree, Tree_CF, CF_Chains, CF_EpiCurves, CF_EpiCurves_Deaths;
	string MetaData; 
	string InitParamValues;

	void init(string scenarioNameArg, string inputFileNameArg, string DJLscenarioNameArg)
	{
		//// record argument variables. 
		inputFileName			= inputFileNameArg; 
		Simon_scenarioName		= scenarioNameArg;
		DJL_scenarioName		= DJLscenarioNameArg; 

		pathInput	= "Data\\"	;
		pathOutput	= "Output\\";
		inputFile	= pathInput + inputFileName;

		string CombinedScenarioName = Simon_scenarioName + DJL_scenarioName;

		Chains				= pathOutput	+ "Chains"				+ CombinedScenarioName + ".txt";
		IndividualR			= pathOutput	+ "individualR"			+ CombinedScenarioName + ".txt";
		Infector			= pathOutput	+ "infector"			+ CombinedScenarioName + ".txt";
		ClusterR			= pathOutput	+ "clusterR"			+ CombinedScenarioName + ".txt";
		ClusterReductionR	= pathOutput	+ "clusterReductionR"	+ CombinedScenarioName + ".txt";
		SimulCluster		= pathOutput	+ "simulCluster"		+ CombinedScenarioName + ".txt";
		Tree				= pathOutput	+ "Trees"				+ CombinedScenarioName + ".txt";
		Tree_CF				= pathOutput	+ "Trees_CF"			+ CombinedScenarioName + ".txt";
		CF_Chains			= pathOutput	+ "CF_Chains"			+ CombinedScenarioName + ".txt";
		CF_EpiCurves		= pathOutput	+ "CF_EpiCurves"		+ CombinedScenarioName + ".txt";
		CF_EpiCurves_Deaths	= pathOutput	+ "CF_EpiCurves_Deaths"	+ CombinedScenarioName + ".txt";
		MetaData			= pathOutput	+ "MetaData"			+ CombinedScenarioName + ".txt";
		InitParamValues		= pathOutput	+ "InitParamValues"		+ CombinedScenarioName + ".txt";
	}
};
struct AllOutput {

	ofstream Chains					; bool Write_Chains					= false;
	ofstream IndividualR			; bool Write_IndividualR			= false;
	ofstream ClusterR				; bool Write_ClusterR				= false;
	ofstream ClusterReductionR		; bool Write_ClusterReductionR		= false;
	ofstream SimulCluster			; bool Write_SimulCluster			= false;
	ofstream Trees					; bool Write_Trees					= false;
	ofstream Trees_CF				; bool Write_Trees_CF				= true;
	ofstream CF_Chains				; bool Write_CF_Chains				= true;
	ofstream CF_EpiCurves			; bool Write_CF_EpiCurves			= true;
	ofstream CF_EpiCurves_Deaths	; bool Write_CF_EpiCurves_Deaths	= true;
	ofstream MetaData				; bool Write_MetaData				= true;
	ofstream InitParamValues		; bool Write_InitParamValues		= true;
		
	void init(FileStrings_Struct FileStrings)
	{
		if (Write_Chains				) Chains				.open(FileStrings.Chains				);
		if (Write_IndividualR			) IndividualR			.open(FileStrings.IndividualR			);
		if (Write_ClusterR				) ClusterR				.open(FileStrings.ClusterR				);
		if (Write_ClusterReductionR		) ClusterReductionR		.open(FileStrings.ClusterReductionR		);
		if (Write_SimulCluster			) SimulCluster			.open(FileStrings.SimulCluster			);
		if (Write_Trees					) Trees					.open(FileStrings.Tree					);
		if (Write_Trees_CF				) Trees_CF				.open(FileStrings.Tree_CF				);
		if (Write_CF_Chains				) CF_Chains				.open(FileStrings.CF_Chains				);
		if (Write_CF_EpiCurves			) CF_EpiCurves			.open(FileStrings.CF_EpiCurves			);
		if (Write_CF_EpiCurves_Deaths	) CF_EpiCurves_Deaths	.open(FileStrings.CF_EpiCurves_Deaths	);
		if (Write_MetaData				) MetaData				.open(FileStrings.MetaData				);
		if (Write_InitParamValues		) InitParamValues		.open(FileStrings.InitParamValues		);
	}
};

enum class VaccCampaignStrategy { 
	REACTIVE	/*Vaccination occurs after first case in hospital, region or national (see ReactiveLevel), plus implementation delay*/, 
	PROACTIVE	/*Vaccination occurs some time before outbreak starts.*/ 
}; // if changing, change Convert_VaccCampaignStrategy_FromString function too. 

enum class ReactiveLevel {
	HOSPITAL, // vaccination reactive to hospital outbreak
	REGIONAL, // vaccination reactive to regional outbreak
	NATIONAL  // vaccination reactive to national outbreak
}; // if changing, change Convert_ReactiveLevel_FromString function too. 

struct ModelRun { //// Set of housekeeping variables

	bool UseCommandLine		= false; // Reading in paramter file from the command line (i.e. with UseCommandLine == true) will overide the parameters below. Otherwise can set them here. 
	bool DJL_InputData		= true;

	int NumIterations			= 11000; 
	int BurnIn					= 1000; 
	int NumCFQuantities			= 0;	//// add to this as required
	int NumCFEpiCurvestStats	= 4;	//// number of summary statistics for counterfactual epidemic curves (mean, median, lower CrI, upper CrI)
	int OutputTreesEvery		= 100;	//// also output counterfactual trees every
	int NumCFsPerTree			= 1;	//// Number counterfactual trees per actual tree. 
	int StoreEvery				= 5;
	std::vector<vector<double>> CF_Chains;				//// 1) iteration; 2) Quantity; vector of vector of doubles. 
	std::vector<vector<double>> CF_EpiCurves;			//// 1) week; 2) Quantity; vector of vector of doubles (need double for mean). 
	std::vector<vector<double>> CF_EpiCurves_Deaths;	//// 1) week; 2) Quantity; vector of vector of doubles (need double for mean). 
	int ** CF_EpiCurves_Internal;						//// indexed by i) iteration; ii) week. Internal version of CF_EpiCurves above. That quantity is outputted. This quantity is used during runtime.
	int ** CF_EpiCurves_Deaths_Internal;				//// indexed by i) iteration; ii) week. Internal version of CF_EpiCurves_Deaths above. That quantity is outputted. This quantity is used during runtime.

	std::vector<string> CF_Names;				//// names of counterfactual quantites

	//// counterfactual parameters. Set as inputs and unchanged throughout MCMC runtime. 
	VaccCampaignStrategy VacCampStrategy	= VaccCampaignStrategy::PROACTIVE;
	ReactiveLevel ReactLevel				= ReactiveLevel::HOSPITAL;
	double Efficacy_Current					= 0.5;		 
	double Efficacy_Start					= 0.0;		 
	double VaccineDuration					= 5000	;	//// Current efficacy = Efficacy * exp (-(TimeSinceVaccination) / VaccineDuration). VaccineDuration == 0 => no waning. 
	double TimeSinceVaccination				= 1;		//// In years, relative to Day_0 (date of first case), can be negative if considering proactive vaccination before outbreak. 
	double Coverage							= 1.0;		
	int ImplementationDelay					= 0;		//// measured in days. 
	int ImmunityDelay						= 0;		//// measured in days. 
	
	///// Add in delays corresponding to incubation period
	bool VaccinateAllHumans				= false;	// this will trump any other booleans (except camels!)
	bool Vaccinate_HCW					= true; 
	double Efficacy_CamelControls		= 0.6;		// efficacy of camel control measures (deliberately vague term to encompass possible vaccination, changes in policy or other control measures aimed at limiting contribution from animal reservoir).

	QuantIndices_Struct QIs;		//// index of quantities 
	AllOutput OUPTUT; 
	FileStrings_Struct FileStrings; 

	int max_threads = 1;
	long seed1 = 547838717ul, seed2 = 943517526ul;

	void AddToCF_Quantities(std::vector<string>& Names, string NameOfThisQuantity, int& CumulativeNumQuantities, int& IndexToSet)
	{
		Names.push_back(NameOfThisQuantity);
		IndexToSet = CumulativeNumQuantities++;
	}
	void init()
	{
		std::cout << "Initialize ModelRun structure " << std::endl;

		AddToCF_Quantities(CF_Names, "PropCasesAverted"			, NumCFQuantities, QIs.PropAverted				);
		AddToCF_Quantities(CF_Names, "PropCasesAverted_HCW"		, NumCFQuantities, QIs.PropAverted_HCW			);
		AddToCF_Quantities(CF_Names, "PropCasesAverted_nHCW"	, NumCFQuantities, QIs.PropAverted_nHCW			);
		AddToCF_Quantities(CF_Names, "DeathsAverted_CF"			, NumCFQuantities, QIs.DeathsAverted_CF			);
		AddToCF_Quantities(CF_Names, "DeathsAverted_CF_HCW"		, NumCFQuantities, QIs.DeathsAverted_CF_HCW		);
		AddToCF_Quantities(CF_Names, "DeathsAverted_CF_nHCW"	, NumCFQuantities, QIs.DeathsAverted_CF_nHCW	);
		AddToCF_Quantities(CF_Names, "FinalCaseDate_CF", NumCFQuantities, QIs.FinalCaseDate_CF);
		
		std::cout << "Model Run NumCFQuantities " << NumCFQuantities << endl; 

		// allocate memory for counterfactual chains and epi curves
		CF_Chains			= vector<vector<double>>(NumIterations			, vector<double>(NumCFQuantities		, 0));
		CF_EpiCurves		= vector<vector<double>>((MAX_ONSET_DAY / 7) + 1, vector<double>(NumCFEpiCurvestStats	, 0));
		CF_EpiCurves_Deaths	= vector<vector<double>>((MAX_ONSET_DAY / 7) + 1, vector<double>(NumCFEpiCurvestStats	, 0));

		Allocate_2D_Array(CF_EpiCurves_Internal			, NumIterations, (MAX_ONSET_DAY / 7) + 1);
		Allocate_2D_Array(CF_EpiCurves_Deaths_Internal	, NumIterations, (MAX_ONSET_DAY / 7) + 1);

		fflush(stderr);	fflush(stdout);
	}
};


