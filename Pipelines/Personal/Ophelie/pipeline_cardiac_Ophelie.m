function pipeline_cardiac_Ophelie(subj_path,subj_file)

%% PIPELINE to analyze cardiac ERP and TF on SEEG
% Given a patient folder, it reads all the parameters from a set_cnfg.m
% file and analyzes the SEEG data, including:
% 
% Syntax:  
%    pipeline_cardiac_Ophelie(subj_path,subj_file);
%
% Inputs:
%   subj_path: 
%       pathway to the main folder of the patient.
%
%   subj_file:
%       [Optional] Name of the file with the configuration parameters of
%       the subject. By default it is "cnfg_patient_cardiac"
%
% Outputs:
%   Different files are saved during the pipeline
%
% See also: cnfg_patient_cardiac cnfg_general_cardiac

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Nov. 2020; Last revision: 26-Mar-2026

% Change log:


%% INITIALIZATION

set(0, 'DefaultTextInterpreter','none')

cnfg_general_cardiac

if nargin == 1, subj_file = 'cnfg_patient_cardiac'; end
cd(subj_path) % Go to patient's folder
eval(subj_file) % Load configuration file

% Define outpath
if ~isfield(cnfg,'outpath')
    cnfg.outpath = [pwd '\Results_cardiac\'];
end
if ~strcmp(cnfg.outpath(end),'\')
    cnfg.outpath = [cnfg.outpath '\'];
end
if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
if ~isfield(cnfg,'infosave'), cnfg.infosave = ''; end
save([cnfg.outpath 'Configuration_file' cnfg.infosave],'cnfg')


%% 1- Load SEEG

% Load SEEG data: 
cfg=[];
cfg.dataset = cnfg.dataset;
%cfg.channel = cnfg.channel;
cfg.channel = {'all', '-ECG'};
cfg.bpfilter = 'yes'; 
cfg.bpfilttype = 'firws';
cfg.bsfilter = 'yes'; 
cfg.bpfreq = cnfg.bpfreq;
cfg.bsfreq = [50 100];
ftdata_cont = ft_preprocessing(cfg);
ftdata_cont.label=correct_elec_names(ftdata_cont.label);
% Create trials:
cfg=[];
cfg.dataset = cnfg.dataset;
cfg.trialdef.eventvalue = cnfg.eventvalue_cardiac;
cfg.trialdef.eventtype = cnfg.eventtype_cardiac; 
cfg.trialdef.prestim    = cnfg.prestim; %in seconds
cfg.trialdef.poststim   = cnfg.poststim; 
cfg = ft_definetrial(cfg);
cfg_trial = cfg;
ftdata_cardiac = ft_redefinetrial(cfg,ftdata_cont);


%% 2- Baseline correction
% Substract baseline for each trial
% All negative values are considered baseline (We can change this and add a
% parameter).

t1 = find(ftdata_cardiac.time{1}<0,1);
t2 = find(ftdata_cardiac.time{1}>0,1);

for trl = 1:size(ftdata_cardiac.trial,2)
    ftdata_cardiac.trial{trl} = ftdata_cardiac.trial{trl} - mean(ftdata_cardiac.trial{trl}(:,t1:t2),2);
end

%% 4- Plot ERPs - Not active. Use only if there are specific channels to check

if false
    relevant_ch = 1:length(ftdata_cardiac.label);

    cfg            = [];
    cfg.channel    = relevant_ch;
    cfg.dosave     = 1;
    cfg.outpath    = [cnfg.outpath 'Plot_ERP\'];
    cfg.plotfig    = 1;
    cfg.trigger    = unique(ftdata_cardiac.trialinfo);
    %cfg.infosave   = cnfg.infosave_cond1;
    cfg.ampyaxis   = 'individual';
    plot_erp(cfg,ftdata_cardiac)
end

%% 5a- COMPONENTS WITH A SIGNIFICANT RESPONSE - NORMAL

if ~isfield(cnfg,'signtrig'), cnfg.signtrig=false; end
if cnfg.signtrig
    disp('Analyzing significant components...')
    
    cfg=[];
    cfg.minsamples = cnfg.signtrig_cfg.minsamples;
    cfg.latency    = cnfg.signtrig_cfg.latency;
    cfg.dosave     = cnfg.dosave;
    cfg.outpath    = [cnfg.outpath 'Evok_resp\'];
    cfg.plotfig    = cnfg.plotfig;
    cfg.infosave   = cnfg.infosave;
    evok_resp = find_sign_trig(cfg,ftdata_cardiac);

    % Select all channels that respond at least to one trigger
    sig_ch = sort(unique([evok_resp.ICsig{:}]));
end

%% 5b- COMPONENTS WITH A SIGNIFICANT RESPONSE - SURROGATES

if ~isfield(cnfg,'signtrig'), cnfg.signtrig=false; end
if cnfg.signtrig
    disp('Analyzing significant components with surrogates...')
    
    cfg=[];
    cfg.Nsurro     = 100;
    cfg.cfg_trial  = cfg_trial;
    cfg.minsamples = cnfg.signtrig_cfg.minsamples;
    cfg.latency    = cnfg.signtrig_cfg.latency;
    cfg.dosave     = cnfg.dosave;
    cfg.outpath    = [cnfg.outpath 'Evok_resp\'];
    cfg.plotfig    = cnfg.plotfig;
    cfg.infosave   = cnfg.infosave;
    cfg.plot_lfdr  = cnfg.plot_lfdr;
    erp_surro(cfg,ftdata_cont);

    % Select all channels that respond at least to one trigger
    sig_ch = sort(unique([evok_resp.ICsig{:}]));
end


%% 7-ITPC and others

if false
if cnfg.tf
    trig = unique(ftdata_cardiac.trialinfo);
    relevant_ch = sig_ch; 
    for ev=1:length(trig)
        cfg=[];
        cfg.channel     = relevant_ch;
        cfg.metric      = {'TF','ITPC'};
        cfg.toi         = -cnfg.prestim_pt:cnfg.tf_cfg.t_step:cnfg.poststim_pt;
        cfg.pval        = cnfg.tf_cfg.pval;
        cfg.Nsurro      = cnfg.tf_cfg.Nsurro;
        cfg.foilim      = cnfg.tf_cfg.foilim;
        cfg.trigger     = trig(ev);
        cfg.doplot      = 0;
        cfg.dosave      = 1;
        cfg.baseline    = [-cnfg.prestim_pt 0];
        cfg.outpath     = [cnfg.outpath 'Time_Frequency\'];
        cfg.infosave    = [cnfg.infosave_cond1 '_trial' num2str(trig(ev))];
        %cfg.trl         = cnfg.trl_pt;
        %cfg.timetrl     = ftdataIC_pt_all.time{1};
        ITPC_MEG_SEEG(cfg,ftdata_pt_tf);
    end
end
end










