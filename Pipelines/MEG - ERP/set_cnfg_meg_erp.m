
%% SET PARAMETERS FOR SUBJECT: SCALES_06  

cnfg = [];

%%% DATA MEG in brainvision (.vhdr)
cnfg.dataset     = 'Data_test\old_new_treated.vhdr';
cnfg.eventtype   = 'Stimulus';
cnfg.eventvalue  = [{'S 520'} {'S 528'} {'S 544'}];
cnfg.chtype      = 'all';

%%% ICA
% ICname is optional. You can compute it the first time and then copy the IC file here:
%cnfg.ICname      = 'D:\Subjects\AScales\Patients\scales_06\MEG\scales_30\01%23%20@02_12\old_new\old_new_treated_sobi_1Hz_90Hz_100.mat';

% If you want to compute ICA:
cnfg.computeICA = true;
cnfg.ica_method = 'sobi'; % 'runica', 'fastica', 'sobi'
cnfg.Ncomp      = 50; % 'all', Number

%NOTE: if computeICA is true and ICname is selected, a new ICA WON'T be computed

%%% PARAMETERS
cnfg.outpath       = 'D:\Subjects\AScales\Patients\scales_06\Treatment\VLM\TEST'; %The pipeline will create subfolders
cnfg.prestim       = 0.5; 
cnfg.poststim      = 1.5;
cnfg.bpfreq        = [1 45];
cnfg.plotfig       = true;
cnfg.dosave        = true;

%% ANALYSES
% Each analysis can be on/off if you want to do it or not

cnfg.signtrig = true; %Significant triggers
    cnfg.signtrig_cfg.minsamples = 10; %Number of consecutive sign. samples
    cnfg.signtrig_cfg.latency    = [0 1]; %Time window to analyze in s
  
cnfg.signcmp  = true; %Compare triggers
    cnfg.signcmp_cfg             = [];
    cnfg.signcmp_cfg.stats       = 'lfdr';
    cnfg.signcmp_cfg.minsamples  = 2; % Only if stats='lfdr';
    cnfg.signcmp_cfg.latency     = [0 1]; %Time window to analyze
    cnfg.signcmp_cfg.trigcmp     = [520 528; 520 544; 528 544]; %triggers to compare 
    
cnfg.ploterp  = false; %Plot selected ERPs, independently of the significance    
    cnfg.ploterp_cfg.channel     = [1 2 3];
    cnfg.ploterp_cfg.trigger     = [528 544];
    cnfg.ploterp_cfg.latency     = [-0.1 1];
    
cnfg.tf = false; %Time Frequency analysis
    cnfg.tf_cfg.channel     = [1 5]; %Channels to analyze. 
    cnfg.tf_cfg.toi         = 0:0.01:1; %Time_start : step : time_end
    cnfg.tf_cfg.foilim      = [1 45]; %Frequencies to analyze
    %[Optional: to correct time-frequency map using the baseline]
    cnfg.tf_cfg.baseline_correction = [-0.3 0]; %[Baseline_start baseline_end]
        
cnfg.CFC = false;     
    cnfg.CFC_cfg.channel = [1 5]; %Channels to analyze.
    cnfg.CFC_cfg.latency = [0 1];
    cnfg.CFC_cfg.f_phase = [3 16 1 2]; %[f_min, f_max, f_step, BW] 
    cnfg.CFC_cfg.f_amp   = [20 120 5 20]; %[f_min, f_max, f_step, BW] 
    cnfg.CFC_cfg.Nsurro  = 10; %Recomended value >100   
    %cnfg.CFC_cfg.ch_interCFC = {'ICA 1' 'B__4'; 'ICA 2' 'B__1'};
    %cnfg.CFC_cfg.ch_interCFC = [1 4; 2 5];
    
cnfg.PS = false;     
    cnfg.PS_cfg.latency = [0 1];    
    
    