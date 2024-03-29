# MERS-CoV vaccine impact counterfactual tree model

This repository contains all code for the Bayesian counterfactual analysis of 
potential MERS-CoV vaccine impact, 
developed by Daniel J. Laydon, and Neil M. Ferguson of the MRC Centre
for Global Infectious Disease Analysis at Imperial College, London, and Simon Cauchemez of the Institut Pasteur, Paris. 
The model assesses vaccine impact by inferring transmission trees, 
("who-infected-whom") analysis, and then "pruning" these trees to generate counterfactuals. 
The repository also contains code to assess and compare various vaccination campaign strategies.
Full details are available at <https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(23)00117-1/fulltext>.


## Building

The model is written in C++. 
Source and header files are located in the subdirectory [MERS_VacTrees](./MERS_VacTrees), 
ready for loading and building from Visual Studio. R output processing code can be found in the 
subdirectory [R](./R), which also contains scripts to build batch files of large cluster jobs. 

## Data

The directory [Data](./Data) contains the anonymised line list data required to replicate our results.
If running elsewhere (e.g. a cluster), then please ensure there is a folder named `Data` within the directory
where the executable is located. 

## Running

The target executable of the C++ source code is named `MERS_VacTrees.exe`, which takes a single 
command line argument specifying a parameter file. The syntax for running the executable is 

`MERS_VacTrees.exe Params_ParticularModelSettings.txt`. 

The code runs reasonably fast on a single core (between approximately 4 and 6 minutes 
for 50,000 iterations when run on the dataset provided). 
However the results of our analysis comprise several thousand 
individual model runs, and so running multiple jobs in parallel on a high-performance cluster is strongly recommended, although
 the model can also be run locally if required. 
If the source code is recompiled having 
set one variable, `UseCommandLine` in `Structs::ModelRun::UseCommandLine` to `false`, then a parameter file is not 
required (this is mostly useful for debugging). Otherwise if `UseCommandLine` is set to `true`, then 
parameter values in `Structs::ModelRun` will be overwritten by values in the parameter file.

Compiled executables will output to the [x64/Release](./x64/Release) folder 
(although this can be changed by specifying a different platform and configuration in Visual Studio, 
or by changing the project properties). The actual executable is not tracked by git.
 
## Parameter files

The executable reads in a parameter file, named `Params_ParticularModelSettings.txt`, 
that govern which features are turned on or off.
Parameter files are generated using the script [MakeBatchAndParamFiles.R](./R/MakeBatchAndParamFiles.R).

The features specified in the parameter files 
produce an output string that identifies the features present in that model run. 
This output string is then included within all output file names. 
The output strings cannot be interpretted without world-class guesswork, or reference to either the C++ function 
`ChooseScenarioName` in [InputOutput.cpp](./MERS_VacTrees/InputOutput.cpp), 
or the R function `ChooseOutputString` in [DirectoryFunctions.R](./R/DirectoryFunctions.R).

We have provided example batch (`.bat`) files to reproduce "proactive" scenarios (`Proactive_WithCamels_Many.bat`)
and "reactive" scenarios (`Reactive_WithCamels_Many.bat`) in the directory [ParamFiles](./ParamFiles), and others can 
be made in the script [MakeBatchAndParamFiles.R](./R/MakeBatchAndParamFiles.R)

## Output

All model outputs are `.txt` files. 
If run locally, model outputs are by default stored in 
[Output](./Output), although outputs are untracked by git as collectively
they can be quite large. If running elsewhere (e.g. a cluster), then ensure there is a folder named `Output` within the directory
where the executable is located.

The folder [R](./R) contains scripts to process and plot model output. By default, these scripts
expect model output to be located in [Output](./Output). 
They further use the R package `here` and expect the working
directory to be the top level of the git repository 
(i.e. whereever the .git is located and the top-level of the Visual studio project).
Plots will be stored in [Plots](./Plots), although these are again untracked by git.

## Workflow

- Compile executable (`MERS_VacTrees.exe`).
- Create batch file and parameter files for collection of jobs/model runs using [MakeBatchAndParamFiles.R](./R/MakeBatchAndParamFiles.R).
- Ensure above batch file is located in same folder as executable. If it isn't, copy it there. Run batch file.
- When finished, summarise model output using [MakeOutputSummaryTable.R](./R/MakeOutputSummaryTable.R).
- If desired, plot posteriors and chains of counterfactuals for individual model runs using [MakeIndividualRunPlots.R](./RMakeIndividualRunPlots.R).
- Plot results using R scripts with `Fig_`, `Figs_` or `Plot` prefixes.

### Relevant papers

The following papers are relevant to the model. Please note that some of them
may require a subscription.

- <https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(23)00117-1/fulltext>
- <https://www.pnas.org/content/113/32/9081>
- <https://www.science.org/doi/10.1126/sciadv.aba8399>

## Copyright and Licensing

Source code is licensed under the GPLv3.

It is Copyright Imperial College of Science, Technology and Medicine. The
lead developers are Daniel J. Laydon, Simon Cauchemez and Neil M. Ferguson. 

This repository includes code modified from
[RANLIB](https://people.sc.fsu.edu/~jburkardt/c_src/ranlib/ranlib.html) which
is licensed under the LGPLv3.

