function pipeline_ICAonSEEG(subj_path,subj_file)

%% PIPELINE to analyze ICA on SEEG (auditory)
% Given a patient folder, it reads all the parameters from a set_cnfg.m
% file and analyzes the SEEG data, including:
%   1- Computes ICA on SEEG data
%   2- Triggers with a significant ERP
%   3- Comparison between triggers
%   4- ITPC
%   5- Plot results
% 
% Syntax:  
%    pipeline_ICAonSEEG(subj_path,subj_file);
%
% Inputs:
%   subj_path: 
%       pathway to the main folder of the patient.
%
%   subj_file:
%       [Optional] Name of the file with the configuration parameters of
%       the subject. By default it is "set_cnfg_ICAonSEEG"
%
% Outputs:
%   Different files are saved during the pipeline
%
% See also: set_cnfg_ICAonSEEG cnfg_general_ICAonSEEG

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Nov. 2020; Last revision: 16-Jan-2026

% Change log:
% 2026-01-23: Added automatic plot of all ERPs (not just significant
% 2026-01-16: Changed the used of loadSEEGdata for standard FieldTrip
% functions

%% INITIALIZATION

set(0, 'DefaultTextInterpreter','none')

cnfg_general_ICAonSEEG

if nargin == 1, subj_file = 'set_cnfg_ICAonSEEG'; end
cd(subj_path) % Go to patient's folder
eval(subj_file) % Load configuration file

% Define outpath
if ~isfield(cnfg,'outpath')
    cnfg.outpath = [pwd '\Results_ICAonSEEG\'];
end
if ~strcmp(cnfg.outpath(end),'\')
    cnfg.outpath = [cnfg.outpath '\'];
end
if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
if ~isfield(cnfg,'infosave'), cnfg.infosave = ''; end
save([cnfg.outpath 'Configuration_file' cnfg.infosave],'cnfg')


%% 1- Load ICA file or compute ICA
%Check for ICA
if isfield(cnfg,'ICname') && cnfg.computeICA
    warning('There is a previous ICA file selected, but a new ICA has been requested. ')
    disp('Please, select only one option:')
    disp('1- Previous ICA file')
    disp('2- Compute new ICA')
    disp('0- Cancel pipeline')
    keypressed = getkey;
    c=1;
    while keypressed ~= 48 && keypressed ~= 49 && keypressed ~= 50  
        warning('Please, the instructions are clear')
        disp('1- Previous ICA file')
        disp('2- Compute new ICA')
        disp('0- Cancel pipeline')
        keypressed = getkey;
        c=c+1;
        if c==4
            error('DEFINITELY, THIS IS TOO DIFFICULT FOR YOU')
        end
    end
    if keypressed == 48
        return
    elseif keypressed == 49
        cnfg.computeICA = false;
    elseif keypressed == 50
        cnfg=rmfield(cnfg,'ICname');
    end
end

%%% If there is a previous ICA file, load it:
if exist([cnfg.outpath 'ICA' cnfg.infosave '.mat'],'file')
    cnfg.computeICA = 0;
    cnfg.ICname = [cnfg.outpath 'ICA' cnfg.infosave '.mat'];
end

% Load SEEG data: 
% CONTINUOUS DATA
%[ftdataSEEG,src_ica]=loadSEEGdata(cnfg); % Load SEEG data
%src_ica.time = []; %Save memory
%src_ica.trial = [];
cfg=[];
cfg.dataset = cnfg.datasetSEEG;
cfg.channel = cnfg.channel;
cfg.bpfilter = 'yes'; 
cfg.bpfilttype = 'firws';
cfg.bsfilter = 'yes'; 
cfg.bpfreq = cnfg.bpfreq;
cfg.bsfreq = [50 100];
ftdataSEEG_cont = ft_preprocessing(cfg);
ftdataSEEG_cont.label=correct_elec_names(ftdataSEEG_cont.label);

% TRIALS
cfg=[];
cfg.dataset = cnfg.datasetSEEG;
cfg.trialdef.eventvalue = cnfg.eventvalueSEEG;
cfg.trialdef.eventtype = cnfg.eventtypeSEEG; 
cfg.trialdef.prestim    = cnfg.prestim; %in seconds
cfg.trialdef.poststim   = cnfg.poststim; 
cfg = ft_definetrial(cfg);
ftdataSEEG = ft_redefinetrial(cfg,ftdataSEEG_cont);

% Compute ICA
if cnfg.computeICA
    cfg = [];
    cfg.method = 'runica';
    src_ica = ft_componentanalysis(cfg, ftdataSEEG_cont);
    src_ica = ica_exp_var(src_ica,ftdataSEEG_cont.trial{1});
end

if isfield(cnfg,'ICname')
    cfg=[];
    cfg.ICname = cnfg.ICname;
    ftdataIC=data2ICA(cfg,ftdataSEEG);
    load([cnfg.outpath 'ICA' cnfg.infosave '.mat'])
else
    cfg=[];
    cfg.src_ica = src_ica;
    ftdataIC=data2ICA(cfg,ftdataSEEG);
    
    if cnfg.dosave
        save([cnfg.outpath 'ICA' cnfg.infosave],'src_ica')
    end
end  

%% 2- Plot ERPs

relevant_ch = 1:length(ftdataIC.label);

cfg            = [];
cfg.channel    = relevant_ch;
cfg.dosave     = 1;
cfg.outpath    = [cnfg.outpath 'Plot_ERP\'];
cfg.plotfig    = 1;
cfg.trigger    = unique(ftdataIC.trialinfo);
%cfg.infosave   = cnfg.infosave_cond1;
plot_erp(cfg,ftdataIC)
    
%% 2- COMPONENTS WITH A SIGNIFICANT RESPONSE
if ~isfield(cnfg,'signtrig'), cnfg.signtrig=false; end
if cnfg.signtrig
    disp('Analyzing significant components...')
    
    cfg=[];
    cfg.minsamples = cnfg.signtrig_cfg.minsamples;
    cfg.latency    = cnfg.signtrig_cfg.latency;
    cfg.dosave     = cnfg.dosave;
    cfg.outpath    = [cnfg.outpath 'Evok_resp_SEEG\'];
    cfg.plotfig    = cnfg.plotfig;
    cfg.infosave   = cnfg.infosave;
    evok_resp = find_sign_trig(cfg,ftdataIC);

    % Select all channels that respond at least to one trigger
    sigIC = sort(unique([evok_resp.ICsig{:}]));
    % Select channels that respond to all triggers
    if length(evok_resp.ICsig)==1, sigIC_all = sigIC;
    else
        sigIC_all = intersect(evok_resp.ICsig{1},evok_resp.ICsig{2});
        for i=3:length(evok_resp.ICsig)
        sigIC_all = intersect(sigIC_all,evok_resp.ICsig{2}); end
    end
    disp(['Components responding to any all trigger: ' num2str(sigIC_all)])
end
       
%% 3- DIFFERENCES BETWEEN TRIGGERS
if ~isfield(cnfg,'signcmp'), cnfg.signcmp=false; end
if cnfg.signcmp
    disp('Analyzing differences between triggers.')

    cfg            = [];
    cfg.channel    = sigIC_all;
    cfg.stats      = cnfg.signcmp_cfg.stats;
    cfg.minsamples = cnfg.signcmp_cfg.minsamples;
    cfg.latency    = cnfg.signcmp_cfg.latency;
    cfg.dosave     = cnfg.dosave;
    cfg.infosave   = cnfg.infosave;
    cfg.outpath    = [cnfg.outpath 'Evok_resp_comp_SEEG\'];
    cfg.plotfig    = cnfg.plotfig;
    ICsig_cmp = [];
    for i = 1:size(cnfg.signcmp_cfg.trigcmp,1)
        cfg.eventvalue = cnfg.signcmp_cfg.trigcmp(i,:);
        evok_resp_cmp=cmp_trig(cfg,ftdataIC);
        ICsig_cmp = [ICsig_cmp evok_resp_cmp.channel];
    end
    ICsig_cmp = unique(ICsig_cmp);
    disp(['Components with differences between conditions: ' num2str(ICsig_cmp)])
end

%% 4-ITPC and others

if ~isfield(cnfg,'doitpc'), cnfg.doitpc=false; end
if cnfg.doitpc
    eventvalue = unique(ftdataIC.trialinfo);
    for ev=1:length(eventvalue)
    cfg=[];
    cfg.channel     = ICsig_cmp;
    cfg.metric      = {'ITPC', 'ITLC'};
    cfg.toi         = cnfg.itpc.toi;
    cfg.foilim      = cnfg.itpc.foilim;
    cfg.trigger     = eventvalue(ev);
    %cfg.Nsurro      = 5;
    %cfg.FDR_q       = 0.05;
    cfg.doplot      = cnfg.plotfig;
    cfg.dosave      = cnfg.dosave;
    cfg.outpath     = [cnfg.outpath 'Connectivity\ITC\'];
    cfg.infosave    = ['trig' num2str(cfg.trigger) cnfg.infosave];

    ITPC_MEG_SEEG(cfg,ftdataIC);
    end
end

%% 5-Plot comparison ICA-SEEG

%PARAMETERS:
resol_SEEG = 2; %distance between signals in units of standard deviation
resol_IC   = 4;

% ERP for each condition.
cond = unique(ftdataSEEG.trialinfo); 
for c=1:length(cond)
    cfg=[];
    cfg.trials  = find(ftdataSEEG.trialinfo==cond(c));
    ftdata_aux  = ft_selectdata(cfg, ftdataSEEG);
    cfg=[];
    cfg.keeptrials = 'yes';
    timelock_SEEG{c} = ft_timelockanalysis(cfg, ftdata_aux);

    cfg=[];
    cfg.trials  = find(ftdataIC.trialinfo==cond(c));
    ftdata_aux  = ft_selectdata(cfg, ftdataIC);
    cfg=[];
    cfg.keeptrials = 'yes';
    timelock_IC{c} = ft_timelockanalysis(cfg, ftdata_aux);
end

% We can put here some conditions to the selected ICAs here... 
cnfg.numICA = min([cnfg.numICA length(sigIC_all)]);
ICAtoplot = sort(sigIC_all); 
ICAtoplot = ICAtoplot(1:cnfg.numICA);

% Divide the figure into subplots
F=figure; hold on

% Obtain position to put the legend
h=subplot(2,length(cond)+1,length(cond)+2);
pos = h.Position;
delete(h)

% a- ICA loadings
subplot(2,length(cond)+1,1), hold on
for i=ICAtoplot
    sig = src_ica.topo(:,i)/max(abs(src_ica.topo(:,i)));
    exv = round(src_ica.exp_var(i)*10000)/100;
    dispname = [ftdataIC.label{i} ' - var=' num2str(exv) '%'];
    plot(sig,1:length(src_ica.label),'LineWidth',1.5,'DisplayName',dispname)
end
legend('Position',pos)
plot([0 0],[0 length(src_ica.label)+1],'k')
axis([-1 1 0 length(src_ica.label)+1])
xlabel('Weight')
ylabel('Channels')
xticks([-1 0 1])
yticks(1:length(src_ica.label))
yticklabels(ftdataSEEG.label)
title('ICA topography')

% b- SEEG traces
for c=1:length(cond)
    subplot(2,length(cond)+1,c+1), hold on
    sig = squeeze(mean(timelock_SEEG{c}.trial,1));
    var_SEEG = std(sig(:));
    t = timelock_SEEG{c}.time;
    for s=1:length(src_ica.label)
        plot(t,sig(s,:)+resol_SEEG*var_SEEG*s,'k')
    end

    axis([t(1) t(end) 0 (s+1)*resol_SEEG*var_SEEG])
    xlabel('Time (s)')
    yticks((1:length(src_ica.label))*resol_SEEG*var_SEEG)
    yticklabels(ftdataSEEG.label)
    title(['Condition: ' num2str(cond(c))])
end    

% c- ICA traces
for c=1:length(cond)
    subplot(2,length(cond)+1,length(cond)+c+2), hold on
    sig = squeeze(mean(timelock_IC{c}.trial(:,ICAtoplot,:),1));
    var_IC = std(sig(:));
    t = timelock_IC{c}.time;
    for i=1:length(ICAtoplot)
        plot(t,sig(i,:)+resol_IC*var_IC*i,'Linewidth',1)
    end

    axis([t(1) t(end) 0 (i+1)*resol_IC*var_IC])
    xlabel('Time (s)')
    yticks((1:length(ICAtoplot))*resol_IC*var_IC)
    yticklabels(ftdataIC.label(ICAtoplot))
end  

if cnfg.dosave
    savefig(F,[cnfg.outpath 'ICAonSEEG' cnfg.infosave])
    saveas(F,[cnfg.outpath 'ICAonSEEG' cnfg.infosave '.png'])
end
if ~cnfg.plotfig
    close(F)
end









