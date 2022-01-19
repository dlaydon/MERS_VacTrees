# MERS-CoV vaccine impact counterfactual tree model

This repository contains all code for the Bayesian counterfactual analysis of the impact of 
a potential MERS-CoV vaccine, developed primarily by Daniel J. Laydon, and Neil M. Ferguson of the MRC Centre
for Global Infectious Disease Analysis at Imperial College, London, and Simon Cauchemez of the Institut Pasteur, Paris. 
The model assesses vaccine impact through the generation and "pruning" of inferred transmission trees, 
("who-infected-whom") analysis. The repository also contains code to compare various vaccination campaign strategies.
Full details are available at ADD LINK.


## IMPORTANT NOTES

:warning: Parameter files and code are named quite esoterically and may be difficult to interpret. 
Updating the code to make it more easily interpretable and user-friendly is a work-in-progress.

## Building

The model is written in C++. 
C++ source and header files are located in the subdirectory [MERS_VacTrees](./MERS_VacTrees), 
ready for loading and building from Visual Studio. R output processing code can be found in the 
subdirectory [R](./R), which also contains scripts to build batch files of large cluster jobs. 

## Data

The directory [Data](./Data) contains the anonymised line list data required to replicate our results. 

## Running

The target executable of the C++ source code is named `MERS_VacTrees.exe`, which takes a single 
command line argument specifying a parameter file. The syntax for running the executable is 

`MERS_VacTrees.exe Params_ParticularModelSettings.txt`. 

The code runs reasonably fast on a single core. However the results of our analysis comprise several thousand 
individual model runs, and so running multiple jobs on a high-performance cluster is strongly recommended. 
However the model can also be run locally. 
If the source code is recompiled having 
set one variable, `UseCommandLine` in `Structs::ModelRun::UseCommandLine` to `false`, then a parameter file is not 
required (this is mostly useful for debugging). Otherwise if `UseCommandLine` is set to `true`, then 
parameter values in `Structs::ModelRun` will be overwritten by values in the parameter file.

For convenience, compiled executables will output to the [x64/Release](./x64/Release) folder 
(although this can be changed by specifying a different platform and configuration in Visual Studio, 
or by changing the project properties)
 
We have provided example batch (`.bat`) files to reproduce an particular "proactive" scenario (`Example_Proactive.bat`)
and an example "reactive" scenario (`Example_Reactive.bat`), 
as well as the script [R/MakeBatchAndParamFiles.R](./R/MakeBatchAndParamFiles.R)

## Parameter files

The executable reads in a parameter file, named `Params_ParticularModelSettings.txt`, 
that govern which features are turned on or off.
Example parameter files are provided in the directory [ParamFiles](./ParamFiles).

By default, the features specified in the parameter files 
produce an output string that identifies the features present in that model run. 
This output string is then included within all output file names. 
The output strings are long, unweildly and difficult to interpret. We apologise for the inconvenience
and are working to make this less painful. 

Our main model (with both age and serotype effects) has output string
`VAC_SILENT_PASSIVE_nENWX_SS_VEs_FSKs3_AGSVEheteroAdd_AS_Hazmult_fAdjHaz_prs1_2_SFU` (called with 
`Example_MainModel.bat`). 
Model with age effects only has output string 
`VAC_SILENT_PASSIVE_nENWX_FSKs3_AGSVEhetero_AS_Hazmult_fAdjHaz_prs1_2_SFU` (called with 
`Example_AgeEffectsOnly.bat`).
Model with serotype effects only has output string 
`VAC_SILENT_PASSIVE_nENWX_SS_VEs_FSKs3_fAdjHaz_prs1_2_SFU` (called with 
`Example_SerotypeEffectsOnly.bat`).
Model with neither age nor serotype effects has output string 
`VAC_SILENT_PASSIVE_nENWX_FSKs3_fAdjHaz_prs1_2_SFU` (called with 
`Example_NoAgeOrSerotypeEffects.bat`).

The above models all assume the "vaccine-as-silent-infection" mechanism, and do not separate cases by disease severity. 
The parameter file, 
`Params_K_SEROPOS_PASSIVE_nENWX_SS_VEs_FSKs3_AGSVEheteroAdd_AS_Hazmult_fAdjHaz_prs1_2_SFU.txt` 
(called with `Example_NoImmunePriming.bat`), models a variant without such immune priming. 


Finally, parameters that model hospitalised and non-hospitalised disease separately are given in 
file `Params_VAC_SILENT_PASSIVE_nENW_MILDSEVERE_hospX_SS_VEs_FSKs3_AGSVEheteroAdd_AS_Hazmult_fAdjHaz_prs1_2_SFU.txt`. 
(called with `Example_Hosp.bat`).
Modelling severe and non-severe disease separately can be accomplished using 
`Params_VAC_SILENT_PASSIVE_nENW_MILDSEVEREX_SS_VEs_FSKs3_AGSVEheteroAdd_AS_Hazmult_fAdjHaz_prs1_2_SFU.txt` 
(called with `Example_Severe.bat`).
Please note that these parameter files use different simulated data files `SimData_Hosp.txt` and `SimData_Severe.txt`
Default output strings are created in the C++ function 
`StructureDefs.h::Housekeeping_Struct::CreateOutputString`, 
and in the R function `DirectoriesEtc.R::ChooseOutputString`.

:warning: All parameter files provided have a relatively small number of MCMC iterations 
(11,000 with 1,000 iteration burn-in) for ease of use. Our real estimates use 1,100,000 iterations with 
a burn-in period of 100,000 iterations. These values can be changed in each individual parameter file
by changing `No_Iterations` and `BurnIn`. 

## Output

All model outputs are `.txt` files. If run locally, model output (e.g. parameter chains, survival tables, attack rates, estimated 
age-specific seroprevalence) is by default stored in 
[ParamFiles/Output](./ParamFiles/Output), although all outputs are untracked by git as 
they can be quite large. 

The folder [R](./R) contains scripts to process and plot model output. By default, these scripts
expect model output to be located in [ParamFiles/Output](./ParamFiles/Output).


### Relevant papers

The following papers are relevant to the model and trial data. Please note that some of them
may require a subscription.

- <https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(12)61428-7/fulltext>
- <https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(14)61060-6/fulltext>
- <https://www.nejm.org/doi/10.1056/NEJMoa1411037?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0www.ncbi.nlm.nih.gov>
- <https://www.nejm.org/doi/full/10.1056/nejmoa1506223>
- <https://science.sciencemag.org/content/353/6303/1033.abstract>

## Copyright and Licensing

Source code is licensed under the GPLv3.

It is Copyright Imperial College of Science, Technology and Medicine. The
lead developers are Daniel J. Laydon and Neil M. Ferguson. 

This repository includes code modified from
[RANLIB](https://people.sc.fsu.edu/~jburkardt/c_src/ranlib/ranlib.html) which
is licensed under the LGPLv3.

