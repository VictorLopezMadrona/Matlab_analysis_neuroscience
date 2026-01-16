%% GENERAL SETTINGS  

%% DATA SEEG
cnfg=[];
cnfg.prestim       = 0.5;
cnfg.poststim      = 1;
cnfg.bpfreq        = [0.5 45];
cnfg.plotfig       = true;
cnfg.dosave        = true;
cnfg.numICA        = 15; %plot up to X ICAs

%% ANALYSES
cnfg.signtrig = true; %Significant triggers
    cnfg.signtrig_cfg.minsamples = 20; %Number of consecutive sig. samples
    cnfg.signtrig_cfg.latency    = [0 1]; %Time window to analyze
    
cnfg.signcmp  = true; %Compare triggers
    cnfg.signcmp_cfg             = [];
    cnfg.signcmp_cfg.stats       = 'lfdr';
    cnfg.signcmp_cfg.minsamples  = 1; %Number of consecutive sig. samples;
    cnfg.signcmp_cfg.latency     = [0 1]; %Time window to analyze
    