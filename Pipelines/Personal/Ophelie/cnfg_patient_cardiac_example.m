%% SET PARAMETERS FOR SUBJECT: SCALES_06  

%This file corresponds to the specific parameters of one subjects.
%It must be in the folder defined in subj_path. If not, the pipeline may
%load a different file.

%% DATA SEEG

cnfg.datasetSEEG    = 'D:\Subjects\AScales\Patients\scales_06\SEEG\3\LogJo0003_treated.vhdr';
%cnfg.ICname
cnfg.eventtypeSEEG  = 'Stimulus';
cnfg.eventvalueSEEG = [{'S 8'} {'S 16'}];

%If not defined, result will be saved on the same location as the cnfg file
%cnfg.outpath        = 'D:\Subjects\AScales\Patients\scales_06\Results_ICAonSEEG\'; %The pipeline will create subfolders

%% ANALYSES    
%%% Compute ICA
cnfg.computeICA  = true; %Compute ICA
cnfg.ica_method  = 'runica'; %runica, fastica, SOBI
cnfg.Ncomp       = 'all';
cnfg.ICApath     = true; %Save ICA data
cnfg.channel     = {'H_*'}; %See FT_CHANNELSELECTION.m for more info

%%% Compare triggers
cnfg.signcmp_cfg.trigcmp     = [8 16];
    
    
    
    