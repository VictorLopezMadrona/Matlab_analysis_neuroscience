
%% Batch file to define all the necessary paths to the toolboxes, parameters and data
% After defining them, the pipeline (flanker_pipeline.m) is computed for
% each subject

%%% NOTE ABOUT FINAL FILE:
% The output is saved for each subject on its own 'results' folder
% It contains 5 values: dERN (all correct-incorrect), dERN (congruency 0%), 
% dERN (congruency 33%), dERN (congruency 66%) and dERN (congruency 100%)

%% DEFINE ALL PATHS HERE
clear 

% Path EEGLAB 
cnfg.eeglab_path = 'C:\Users\lopezmadrona\Desktop\Victor_Lopez\Matlab functions\eeglab2025.0.0';   
% Path FIELDTRIP
cnfg.fieldtrip_path = 'C:\Users\lopezmadrona\Desktop\Victor_Lopez\Matlab functions\fieldtrip-20240110'; 
% Path Matlab_analysis_neuroscience
cnfg.neurotoolbox_path = 'C:\Users\lopezmadrona\Desktop\Victor_Lopez\Matlab functions\Matlab_analysis_neuroscience';

% Path to store the results
% It will create one folder per subject and save the following info: 
% dERN.txt, summary.txt, src_ica_flanker.mat (ICA file for replicability)
cnfg.outpath_all = 'C:\Users\lopezmadrona\Desktop\Victor_Lopez\Projects\EEGMANYANALYSTS\Results';

% Path with all the EEG data
% NOTE: We expect to have here as many folders as participants, and ONLY
% that, i.e., all the folders here should have the format
% \sub-XXX\eeg\sub-xxx_task-Flanker_eeg.set and so on
subj_path = 'C:\Users\lopezmadrona\Desktop\Victor_Lopez\Projects\EEGMANYANALYSTS\Subsample_10';

%% Find all the participants

folders = dir(subj_path);
subj_names = {folders([folders.isdir] & ~startsWith({folders.name}, '.')).name};

% Check that all the folders in this path start with 'sub-', i.e., they
% all belong to different participants
allStartsWithSub = all(startsWith(subj_names, 'sub-'));
if ~allStartsWithSub
    error(['The folders in ' subj_path ' do not have the expected format'])
end

%% Initialization

%%% SOME PARAMETERS

% Trial duration
cnfg.stimdef = [-2 0.5]; 
% Filter
cnfg.filter = [];
% High-pass filter at 0.25 Hz using FIR (more stable for low frequencies)
cnfg.filter.hpfilter = 'yes';
cnfg.filter.hpfreq = 0.25;
cnfg.filter.hpfilttype = 'fir';
cnfg.filter.hpfiltdir = 'twopass';  % Zero-phase filtering for no phase distortion
% Low-pass filter at 45 Hz using FIR
cnfg.filter.lpfilter = 'yes';
cnfg.filter.lpfreq = 45;
cnfg.filter.lpfilttype = 'fir';
cnfg.filter.lpfiltdir = 'twopass';
% Notch (band-stop) filter for 50 and 100 Hz
cnfg.filter.bsfilter = 'yes';
cnfg.filter.bsfreq = [49 51; 99 101];  % Slightly widened stopband
% Detrend data
cnfg.filter.detrend = true;

%%% ADD PATH

% Add EEGLAB to path
addpath(cnfg.eeglab_path);   
% Add FIELDTRIP to path 
addpath(cnfg.fieldtrip_path); 
% Add Matlab_analysis_neuroscience to path
addpath(genpath(cnfg.neurotoolbox_path));

%% Sequentialy, run the pipeline for each participant

% Create results folder if it doesn't exist
if ~exist(cnfg.outpath_all,'dir'), mkdir(cnfg.outpath_all); end

clear dERN_output % THIS CAN BE REMOVED
ec=1; % Counter to store errors
clear encounter_error
for subji = 1:length(subj_names)
    
    % Filename and outpath for the subject
    cnfg.filename = [subj_path '\' subj_names{subji} '\eeg\' subj_names{subji} '_task-Flanker_eeg.set'];
    cnfg.outpath = [cnfg.outpath_all '\' subj_names{subji} '\'];
    
    % Create results folder if it doesn't exists
    if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
    
    % Check if we already have the results for this subject, skip it:
    if ~exist([cnfg.outpath 'dERN.txt'],'file')
        try
            dERN_output(subji,:)=flanker_pipeline(cnfg);
        catch ME    
            encounter_error{ec} = ([subj_names{subji} ' ', ME.message]);
            ec=ec+1;
        end
    end
end
