function [ftdata]=pipeline_meg_erp(cnfg_path)

%% PIPELINE to read and analyze MEG (and SEEG)
%
% NOTE: This pipeline does not preprocess the data. We recommend AnyWave to do so.
%
% Given a patient folder, it reads all the parameters from a set_cnfg_oldnew.m
% file and analyzes the MEG data, including:
%   1- Compute ICA 
%   2- Triggers with a significant ERP
%   3- Comparison between triggers
%   4- Plot selected ERPs
%   5- Time-freq analysis and ITPC 
%   6- Cross-Frequency Coupling
%   7- Power Spectrum
%
% Syntax:  
%    pipeline_meg_erp(cnfg_path);
%
% Inputs:
%   cnfg_path: 
%       pathway and filename of the cnfg file of the patient
%
% Outputs:
%   Different files are saved during the pipeline
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug. 2021; Last revision: 12-Aug-2024

% Change log:
% 12/08/2024: Correct filtering
% 28/10/2021: Compute ICA  
% 13/08/2021: Included time-frequency analysis and ITPC

% TO DO:


%% INITIALIZATION

set(0, 'DefaultTextInterpreter','none')

[path_cfg,name_cfg,~] = fileparts(cnfg_path);
cd(path_cfg)
eval(name_cfg) % Load configuration file

if ~strcmp(cnfg.outpath(end),'\')
    cnfg.outpath = [cnfg.outpath '\'];
end
if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
save([cnfg.outpath 'Configuration_file'],'cnfg')

%% Load data

%% LOAD RAW-DATA

freq_filt = cnfg.bpfreq;

cfg=[];
cfg.dataset = cnfg.dataset;
cfg.channel = cnfg.chtype;
cfg.bpfilter = 'yes'; 
cfg.bpfilttype = 'firws';
cfg.bpfreq = freq_filt;
cfg.bsfilter = 'yes'; 
cfg.bsfreq = [49.5 50.5];
ftdata = ft_preprocessing(cfg);

cfg=[];
cfg.dataset = filename;
cfg.trialdef.eventvalue = cnfg.eventvalue;
cfg.trialdef.eventtype = cnfg.eventtype;
cfg.trialdef.prestim = cnfg.prestim; % in seconds
cfg.trialdef.poststim = cnfg.poststim; % in seconds
ftdata = ft_redefinetrial(cfg,ftdata);


% Remove bad channels marked in anywave
[data_path,data_name,~] = fileparts(cnfg.dataset);
bad_ch_fullfile = [data_path '\' data_name '.bad'];
if exist(bad_ch_fullfile,'file')
    fileID = fopen(bad_ch_fullfile);
    bad_ch_labels = textscan(fileID,'%s','delimiter','\n');
    bad_ch_labels = bad_ch_labels{1};
    fclose(fileID);
    rm_channels = strcat('-', bad_ch_labels)'; %Substract bad channels
    rm_channels = [{'all'}, rm_channels(:)'];
    cfg = [];
    cfg.channel  = rm_channels;
    ftdata = ft_selectdata(cfg, ftdata);
end

% If ICA file, obtain ICA components
if isfield(cnfg,'ICname')
    cnfg.computeICA = false;
    ftdata = data2ICA(cnfg,ftdata); % Load IC data
end

% Remove artefacts selected in AnyWave
[filepath,name,~] = fileparts(cnfg.dataset);
cfg=[];
cfg.filename = [filepath '\' name, '.mrk'];
if exist(cfg.filename,'file')
    ftdata = remove_trl_art(cfg,ftdata);
    cfg.marker_name = {'artifact'};
    ftdata = remove_trl_art(cfg,ftdata);    
else
    cfg.filename = [filepath '\' name, '.vhdr.mrk'];
    ftdata = remove_trl_art(cfg,ftdata);
    cfg.marker_name = {'artifact'};
    ftdata = remove_trl_art(cfg,ftdata);
end

%% 1- COMPUTE ICA

if ~isfield(cnfg,'computeICA'), cnfg.computeICA=true; end
if cnfg.computeICA
    
    if ~isfield(cnfg,'ica_method'), cnfg.ica_method='runica'; end
    if ~isfield(cnfg,'Ncomp'), cnfg.ica_method='all'; end
    
    % 1- Concatenate trials
    %ftdataMEG=ft_ARIMA(ftdataMEG); %Pre-whiten power spectrum
    ft_conc_data = concatenate_fttrials(ftdata);
    if strcmp(cnfg.Ncomp,'all')
        cnfg.Ncomp = length(ftdata.label);
    end
    
    % 1.5- Compute PCA?
    if strcmp(cnfg.ica_method,'sobi')
        data_pca = ft_conc_data;      
        [data_pca.trial{1},eigenvectors]=dyn_pca(ft_conc_data.trial{1},cnfg.Ncomp);
        for n=1:cnfg.Ncomp, label_pca{n}=['PCA_' num2str(n)]; end
        data_pca.label = label_pca';
    else
        data_pca=ft_conc_data;
    end
    
    % 2- Compute ICA
    cfg = [];
    cfg.method       = cnfg.ica_method;
    cfg.numcomponent = cnfg.Ncomp;
    %if isfield(cnfg,'sobi')
    %if isfield(cnfg.sobi,'p_correlations')
    %    cfg.sobi.p_correlations = cnfg.sobi.p_correlations;
    %    cfg.sobi.n_sources = cnfg.Ncomp;
    %end     
    %end
    cfg.demean       = 'no';  % whether to baseline correct the data not demean 
                              % the data as stated  in the doc (default = 'yes')
    src_ica = ft_componentanalysis(cfg, data_pca);

    % 2.5- Undo PCA
    if strcmp(cnfg.ica_method,'sobi')
        src_ica.unmixing  = src_ica.unmixing * eigenvectors(:,1:cnfg.Ncomp)';
        src_ica.topo      = pinv(src_ica.unmixing);
        src_ica.topolabel = ft_conc_data.label;
    end
    
    src_ica = ica_exp_var(src_ica,ft_conc_data.trial{1});
       
    cfg=[];
    cfg.src_ica = src_ica;
    ftdata = data2ICA(cfg,ftdata); % Load IC data
        
    %SAVE ICA
    modality = 'MEG';
    lpf = cnfg.bpfreq(2);
    hpf = cnfg.bpfreq(1);
    sr = ftdata.fsample;
    mixing = src_ica.topo;
    unmixing = src_ica.unmixing;
    labels = src_ica.topolabel;
    exp_var = src_ica.exp_var;
        
    [filepath,name,~] = fileparts(cnfg.dataset);
    outname = [filepath '\' name '_' cnfg.ica_method '_' num2str(hpf) 'Hz_' num2str(lpf) 'Hz_' num2str(cnfg.Ncomp) 'test'];% 'c_trig' name_trig];
    save(outname,'modality','lpf','hpf','sr','mixing','unmixing','labels','exp_var')

end

%% 2- TOPOGRAPHY 

outpath=[cnfg.outpath 'IC_topo\'];
if ~exist(outpath,'dir')
    mkdir(outpath)
end
if isfield(src_ica,'labels') % ICA from fieldtrip or from Anywave
    src_out.fsample = src_ica.sr;
    src_out.time = [];
    src_out.trial = [];
    src_out.topo = src_ica.mixing;
    src_out.unmixing = src_ica.unmixing;
    src_out.label = [];
    src_out.topolabel = src_ica.labels;
    src_out.grad = [];
    src_out.trialinfo = [];
    src_out.cfg = [];
    src_ica = src_out;
end
for nic=1:size(src_ica.topo,2)
    cfg=[];
    cfg.component = nic;
    if ~isfield(cnfg,'layout')
        cfg.layout = '4D248_helmet';
    else
        cfg.layout = cnfg.layout;
    end
    h=figure; ft_topoplotIC(cfg,src_ica); colorbar
    savefig(h,[outpath 'IC_' num2str(nic)])
    saveas(h,[outpath 'IC_' num2str(nic) '.png'])
    close(h);
end


%% 3- COMPONENTS WITH A SIGNIFICANT RESPONSE

if ~isfield(cnfg,'signtrig'), cnfg.signtrig=false; end
if cnfg.signtrig
    disp('Analyzing significant components...')
    
    cfg=[];
    cfg.minsamples = cnfg.signtrig_cfg.minsamples;
    cfg.latency    = cnfg.signtrig_cfg.latency;
    cfg.dosave     = cnfg.dosave;
    cfg.outpath    = [cnfg.outpath 'Evok_resp\'];
    cfg.plotfig   = cnfg.plotfig;
    evok_resp = find_sign_trig(cfg,ftdata);
end

%% 4- DIFFERENCES BETWEEN TRIGGERS

if ~isfield(cnfg,'signcmp'), cnfg.signcmp=false; end
if cnfg.signcmp
    disp('Analyzing differences between triggers...')
    
    cfg            = [];
    cfg.channel    = 'all';
    cfg.stats      = cnfg.signcmp_cfg.stats;
    cfg.minsamples = cnfg.signcmp_cfg.minsamples;
    cfg.latency    = cnfg.signcmp_cfg.latency;
    cfg.dosave     = cnfg.dosave;
    cfg.outpath    = [cnfg.outpath 'Evok_resp_comp\'];
    cfg.plotfig    = cnfg.plotfig;
    for i = 1:size(cnfg.signcmp_cfg.trigcmp,1)
        cfg.eventvalue = cnfg.signcmp_cfg.trigcmp(i,:);
        evok_resp_cmp=cmp_trig(cfg,ftdata);
    end
end
 
%% 5- Plot ERPs

if ~isfield(cnfg,'ploterp'), cnfg.ploterp=false; end
if cnfg.ploterp
    disp('Plotting selected channels...')
    
    cfg            = [];
    cfg.channel    = cnfg.ploterp_cfg.channel;
    cfg.trigger    = cnfg.ploterp_cfg.trigger;
    if isfield(cnfg.ploterp_cfg,'latency')    
        cfg.latency = cnfg.ploterp_cfg.latency; end
    cfg.dosave     = cnfg.dosave;
    cfg.outpath    = [cnfg.outpath 'Plot_ERP\'];
    cfg.plotfig    = cnfg.plotfig;
    plot_erp(cfg,ftdata)
    
end

%% 6- Time-Frequency Analysis (and ITPC)

if ~isfield(cnfg,'tf'), cnfg.tf=false; end
if cnfg.tf
    disp('Analyzing Time-Frequency responses')
    trig = unique(ftdata.trialinfo);
    for ev=1:length(trig)
        cfg=[];
        cfg.channel     = cnfg.tf_cfg.channel;
        cfg.metric      = {'TF','ITPC'};
        cfg.toi         = cnfg.tf_cfg.toi;
        cfg.foilim      = cnfg.tf_cfg.foilim;
        cfg.trigger     = trig(ev);
        cfg.doplot      = cnfg.plotfig;
        cfg.dosave      = cnfg.dosave;
        cfg.baseline_correction = cnfg.tf_cfg.baseline_correction;
        cfg.outpath     = [cnfg.outpath 'Time_Frequency\'];
        cfg.infosave    = ['_trial' num2str(trig(ev))];
        
        ITPC_MEG_SEEG(cfg,ftdata);
    end
end

%% 7- Cross-Frequency Coupling

if ~isfield(cnfg,'CFC'), cnfg.CFC=false; end
if cnfg.CFC
    disp('Analyzing CFC')
    
    trig = unique(ftdata.trialinfo);
    for ev=1:length(trig)
        cfg = [];
        cfg.trials = find(ismember(ftdata.trialinfo,trig(ev)));
        ftdata_aux = ft_selectdata(cfg,ftdata);
        
        cfg = [];
        cfg.channel  = cnfg.CFC_cfg.channel;
        if isfield(cnfg.CFC_cfg,'ch_interCFC')
            cfg.ch_interCFC = cnfg.CFC_cfg.ch_interCFC; end
        cfg.latency  = cnfg.CFC_cfg.latency;
        cfg.f_phase  = cnfg.CFC_cfg.f_phase;
        cfg.f_amp    = cnfg.CFC_cfg.f_amp;
        cfg.Nsurro   = cnfg.CFC_cfg.Nsurro;
        cfg.doplot   = cnfg.plotfig;
        cfg.dosave   = cnfg.dosave;
        cfg.outpath  = [cnfg.outpath 'CFC\']; 
        cfg.infosave = ['_trig' num2str(trig(ev))];
        %comodulogram_ft(cfg,ftdata_aux);
        comodulogram_trial(cfg,ftdata_aux);
    end
end

%% 8- POWER SPECTRUM

if ~isfield(cnfg,'PS'), cnfg.PS=false; end
if cnfg.PS
    disp('Computing PS')
    
    trig = unique(ftdata.trialinfo);
    for ev=1:length(trig)
        cfg = [];
        cfg.latency = cnfg.PS_cfg.latency;
        cfg.trials  = find(ismember(ftdata.trialinfo,trig(ev)));
        ftdata_aux  = ft_selectdata(cfg,ftdata);
        
        cfg = [];
        cfg.doplot   = cnfg.plotfig;
        cfg.dosave   = cnfg.dosave;
        cfg.outpath  = [cnfg.outpath 'PS\']; 
        cfg.infosave = [cnfg.infosave '_trig' num2str(trig(ev))];
        computePS_ft(ftdata_aux,cfg);
    end
end

