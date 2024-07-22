
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ABOUT THE PIPELINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The goal of this pipeline is to analyze MEG experiments based on trials.
It will do an ICA and analyze the responses of each component, including:
ERP, Power Spectrum, time-freq analysis, inter-trial phase coherence (ITPC), Cross-Frequency Coupling...

The pipeline should be run a first time to compute the ICA matrix and identified the components with a significant ERP.
For the analysis based on time-frequency maps, we need to specify the components to analyze, as it is an operation quite slow.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% HOW TO USE THE PIPELINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PREPROCESS DATA

	There is no preprocessing in the pipeline. First, load the data in anywave and preprocess it there (select artefacts, bad channels, filtering...)
	To run the pipeline, the data must be in .vhdr format (brainvision), so save the data in that format.
	

PIPELINE

0- You need to have FieldTrip for Matlab installed (https://www.fieldtriptoolbox.org/)

1- Open matlab

2- You need to add the files and functions to the matlab path:
	2a- In matlab go to HOME -> Set Path -> Add with subfolders
     	2b- Select the main folder, this is, the folder containing all the matlab files.
     	2c- Click on select folder -> save -> close

3- Ok, now you see that there are two matlab files in the folder: "pipeline.m" and "set_cnfg.m"
	The first one is the function to run in matlab.
	The second one contains all the parameters.
	By default, the pipeline should work, but if you want to modify something, that is the file to do it.

4- To run the code just write in the command window:
	>> pipeline_meg_erp(pathway_set_cnfg)
	where pathway_set_cnfg is the pathway and name to the configuration file