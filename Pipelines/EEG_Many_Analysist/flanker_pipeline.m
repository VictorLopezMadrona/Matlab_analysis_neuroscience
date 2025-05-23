function dERN_output = flanker_pipeline(cnfg)

%% PIPELINE to read and analyze data from EEG MANY ANALYSTS

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Mar. 2025; Last revision: 22-May-2025


%% Load continuous data
disp('Loading data...')

% Launch eeglab and close it
eeglab; close

% Load the dataset
EEG = pop_loadset(cnfg.filename);

% Convert to FieldTrip format 
ft_defaults;
ftdata = eeglab2fieldtrip(EEG, 'preprocessing', 'none');

%% Preprocessing ftdata - Filtering
disp('Filtering...')
% Filtering
cfg = [];
cfg = cnfg.filter;
% Apply filters to already loaded data
ftdata = ft_preprocessing(cfg, ftdata);

%% Preprocessing ftdata - Bad channels
disp('Identifying bad channels...')
% Select channels with very high variance and mark them as 'bad'.
% These bad channels influence a lot ICA compared to bad channels with
% low variance.
z_var = zscore(var(ftdata.trial{1}, 0, 2));
bad_var = z_var > 3.5;

BAD_CH = ftdata.label(bad_var); % <--- Info with bad channels

% Select all channels except bad channels
cfg = [];
cfg.channel = setdiff(ftdata.label, ftdata.label(bad_var));  % keep all except bad channels
ftdata = ft_selectdata(cfg, ftdata);

%% Create trials
disp('Creating trials...')

% Select Trials of Interest from EEGLAB
ToI=find(ismember({EEG.event.type},{'107','117','127','137','109','119','129','139'})==1);
latencies = cell2mat({EEG.event.latency});
eventtype = cellfun(@str2double, {EEG.event.type});

% Define trial duration
onset_samples = latencies(ToI);
pre_samples = round(cnfg.stimdef(1)*ftdata.fsample);
post_samples = round(cnfg.stimdef(2)*ftdata.fsample);

% Define trials
trl = zeros(length(ToI), 3);
trl(:, 1) = onset_samples + pre_samples;  % Start sample
trl(:, 2) = onset_samples + post_samples; % End sample
trl(:, 3) = pre_samples;                  % Trial onset after sample starts

% Create ftdata_trial
cfg = [];
cfg.trl = trl;  % Assign the trial matrix
ftdata = ft_redefinetrial(cfg, ftdata);
ftdata.trialinfo = eventtype(ToI)';

%% Remove noisy trials iteratively
disp('Identifying bad trials...')
max_iter = 20;
bad_trials_all = [];
ftdata_aux = ftdata;

% Recursively identify noisy trials and remove them
% The loop allows to recompute zscore after removing trials that are
% outlayers.
for iter = 1:max_iter
    
    % Initialize
    Ntrials = length(ftdata_aux.trial);
    var_vec   = zeros(Ntrials, 1);
    kurt_vec  = zeros(Ntrials, 1);
    
    % Compute variance and kurtosis per trial
    for i = 1:Ntrials
        x = ftdata_aux.trial{i}(:);  % Flatten channels × time to 1D
        var_vec(i)  = var(x);
        kurt_vec(i) = kurtosis(x);
    end
    
    % Z-score both metrics
    z_var  = (var_vec  - mean(var_vec))  / std(var_vec);
    z_kurt = (kurt_vec - mean(kurt_vec)) / std(kurt_vec);
    
    % Set thresholds 
    thresh_var  = 3;  
    thresh_kurt = 3; 
    
    % Flag bad trials
    bad_var  = z_var  > thresh_var;
    bad_kurt = abs(z_kurt) > thresh_kurt;
    
    % Combine criteria
    bad_trials_idx = find(bad_var | bad_kurt);  % union
    
    if sum(bad_trials_idx)==0
        break;
    end
    
    % Store bad trials to keep track
    bad_trials_all = [bad_trials_all bad_trials_idx'];
    
    good_trials_idx = setdiff(1:Ntrials, bad_trials_idx);
    
    % Remove bad trials
    cfg = [];
    cfg.trials = good_trials_idx;
    ftdata_aux = ft_selectdata(cfg, ftdata_aux);
end

ftdata = ftdata_aux;
BAD_TRIAL = length(bad_trials_all); % <--- Number of bad trials

%% COMPUTE OR LOAD ICA

% Add Matlab_analysis_neuroscience to path to be sure that we used the
% correct runica function
addpath(genpath(cnfg.neurotoolbox_path));

% If there is already an ICA file in the outpath folder, load it.
% This is for replicability, as the ICA matrix may be different.
if exist([cnfg.outpath 'src_ica_flanker.mat'],'file')
    disp('Loading ICA matrix from previous file...')
    cnfg.ICname = [cnfg.outpath 'src_ica_flanker.mat'];
end

% If there is no ICA file, compute it using 'runica' and fieldtrip
if ~isfield(cnfg,'ICname')
    disp('Computing ICA...')
    % Fix random seed to obtain same ICA
    rng(123)
    cfg = [];
    cfg.method = 'runica';   
    src_ica = ft_componentanalysis(cfg, ftdata);
    %src_ica = ica_exp_var(src_ica,ftdata.trial{1});
    
    %%% Save ICA file
    src_ica.time = []; %Save space
    src_ica.trial = [];
    save([cnfg.outpath 'src_ica_flanker'],'src_ica')
else
    % Load ICA file if possible
    load(cnfg.ICname,'src_ica')
end

% I use the ICA matrix (computed or loaded) to obtain the components.
% ftdataICA contains the ICA time-series. This is the data that will be
% analyzed in this pipeline
cfg=[];
cfg.src_ica = src_ica;
ftdataICA   = data2ICA(cfg,ftdata);

%% Remove noisy ICA as well
disp('Identifying bad ICs...')
% Select ICs with high variance, peak2peak or kurtosis. Do not analyze
% them.

% Recursively identify noisy trials and remove them
% The loop allows to recompute zscore after removing trials that are
% outlayers.
max_iter = 10;
bad_channels_all = [];
ftdata_aux = ftdataICA;
for iter = 1:max_iter
    conc_data = concatenate_fttrials(ftdata_aux);
    data = conc_data.trial{1};
    
    %Variance
    z_var = zscore(var(data, 0, 2));
    bad_var = abs(z_var) > 3;

    %Peak 2 peak
    p2p = max(data, [], 2) - min(data, [], 2);
    z_p2p = zscore(p2p);
    bad_p2p = abs(z_p2p) > 3;

    %Kurtosis
    z_kurt = zscore(kurtosis(data, 0, 2));
    bad_kurt = abs(z_kurt) > 3;

    bad_mask = bad_var | bad_p2p | bad_kurt;
    
    if sum(bad_mask)==0
        break;
    end
    
    % Store bad channels (if you want to keep track)
    bad_channels_all = [bad_channels_all conc_data.label(bad_mask)'];
    
    % Remove bad channels from ftdata
    cfg = [];
    cfg.channel = setdiff(ftdata_aux.label, ftdata_aux.label(bad_mask));  % keep only good chans
    ftdata_aux = ft_selectdata(cfg, ftdata_aux);
end

[~,bad_channels_idx] = ismember(bad_channels_all,ftdataICA.label);
ftdataICA = ftdata_aux;


%% Z-score each Independent Component separately

% concatenate all trials
conc_data = concatenate_fttrials(ftdataICA);

% Compute mean and std per channel
mu = mean(conc_data.trial{1}, 2);              % [channels × 1]
sigma = std(conc_data.trial{1}, 0, 2);         % [channels × 1]

% Apply z-scoring to each trial
for i = 1:length(ftdataICA.trial)
    ftdataICA.trial{i} = (ftdataICA.trial{i} - mu) ./ sigma;
end

%% Select component with prefrontal topography

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% WHAT IF NO GOOD TOPO? %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only interested in Independent Components with a prefrontal topograhy:
% i.e., the maximal values of the spatial distribution (mixing matrix)
% should be in one of these channels:
ChoI = {'CZ', 'FZ', 'FC1', 'FC2', 'C1', 'C2'};
chm = ismember(ftdata.label, ChoI);
chlist = find(chm==1);

% Obtain spatial topo of all ICAs
ica_topo = src_ica.topo;
% Remove bad ICAs (selected before)
ica_topo(:,bad_channels_idx)=[];

% I select only those component with a prefrontal topography.
% Max in absolute value (i.e., max or min), 
% as the ICA topography doesn't have polarity
ICselect = [];
for chi=1:size(ica_topo,2)
    [~,pmax]=max(ica_topo(:,chi));
    [~,pmin]=min(ica_topo(:,chi));
    
    if ismember(pmax,chlist) || ismember(pmin,chlist)
        % These are the Independent components with 'prefrontal' topography
        ICselect = [ICselect chi];
    end
end

%% Compute ERP difference between conditions

% Create new trials just correct vs incorrect
% I relabel the triggers, naming '7' all the correct trials and '9' all the
% incorrect, independently of the congruency
ftdataERP_aux = ftdataICA;
ftdataERP_aux.trialinfo(ftdataERP_aux.trialinfo==107) = 7;
ftdataERP_aux.trialinfo(ftdataERP_aux.trialinfo==117) = 7;
ftdataERP_aux.trialinfo(ftdataERP_aux.trialinfo==127) = 7;
ftdataERP_aux.trialinfo(ftdataERP_aux.trialinfo==137) = 7;
ftdataERP_aux.trialinfo(ftdataERP_aux.trialinfo==109) = 9;
ftdataERP_aux.trialinfo(ftdataERP_aux.trialinfo==119) = 9;
ftdataERP_aux.trialinfo(ftdataERP_aux.trialinfo==129) = 9;
ftdataERP_aux.trialinfo(ftdataERP_aux.trialinfo==139) = 9;

% Compute ERP in a time-window around the response onset (-100ms to 250ms)
cfg = [];
cfg.channel = ICselect;
cfg.latency = [-0.1 0.25];
cfg.trials = (ftdataERP_aux.trialinfo==7);
tl7 = ft_timelockanalysis(cfg, ftdataERP_aux);
cfg.trials = (ftdataERP_aux.trialinfo==9);
tl9 = ft_timelockanalysis(cfg, ftdataERP_aux);

%% Select best component: ICoI

% I select the component with maximal differences between both conditions
% (correct vs uncorrect) after response onset. I expect the maximal
% response to be between 0 and 100ms after response.
timewindow = [0 0.1]; % Time window of interest in seconds
%Time window in samples:
samplewindow(1) = find(tl7.time>timewindow(1),1);
samplewindow(1) = samplewindow(1)-1; %To include t=0 as well
samplewindow(2) = find(tl7.time>timewindow(2),1);

clear maxval %store here the dERN of each component
clear del_maxval %store here the lag of the max dERN
for chi=1:length(ICselect) %Loop across all component with good topography
    dERN_aux = tl7.avg(chi,:)-tl9.avg(chi,:); %Compute difference of ERP
    %Remove the mean of the difference. This is to focus on differences at
    %the time window of interest.
    dERN_aux = dERN_aux - mean(dERN_aux); 
    %Compute max of the difference in the time-window of interest. Absolute
    %value because I dont know the polarity after ICA
    [maxval(chi),del_maxval(chi)] = max(abs(dERN_aux(samplewindow(1):samplewindow(2))));
end

% Select ICA with maximal differences
[~,p] = max(maxval);
ICoI = ICselect(p); % <---- This is my final selected component

% Select the time-lag of the maximal dERN. 
% dERN will be computed at this time-point
lag_dERN = samplewindow(1) + del_maxval(p) - 1; % <---- LAG


%% Compute dERN on the best component (ICoI)

% Difference of the ERP between two main conditions for selected IC:
dERP = tl7.avg(p,:)-tl9.avg(p,:);
dERN = abs(dERP(lag_dERN)); % <--- THIS IS THE FINAL RESULT

% Compute dERN for each congruency
cfg = [];
cfg.channel = ICoI;
cfg.latency = [-0.1 0.25];

% 0% Congruency
% Control there are trials to compute it
if sum(ftdataICA.trialinfo==107)==0 || sum(ftdataICA.trialinfo==109)==0
    dERN_0 = 0;
else
    cfg.trials = (ftdataICA.trialinfo==107);
    tl7 = ft_timelockanalysis(cfg, ftdataICA);
    cfg.trials = (ftdataICA.trialinfo==109);
    tl9 = ft_timelockanalysis(cfg, ftdataICA);
    dERP = tl7.avg-tl9.avg;
    dERN_0 = abs(dERP(lag_dERN)); % <--- THIS IS THE FINAL RESULT
end

% 33% Congruency
if sum(ftdataICA.trialinfo==117)==0 || sum(ftdataICA.trialinfo==119)==0
    dERN_33 = 0;
else
    cfg.trials = (ftdataICA.trialinfo==117);
    tl7 = ft_timelockanalysis(cfg, ftdataICA);
    cfg.trials = (ftdataICA.trialinfo==119);
    tl9 = ft_timelockanalysis(cfg, ftdataICA);
    dERP = tl7.avg-tl9.avg;
    dERN_33 = abs(dERP(lag_dERN)); % <--- THIS IS THE FINAL RESULT
end

% 66% Congruency
if sum(ftdataICA.trialinfo==127)==0 || sum(ftdataICA.trialinfo==129)==0
    dERN_66 = 0;
else
    cfg.trials = (ftdataICA.trialinfo==127);
    tl7 = ft_timelockanalysis(cfg, ftdataICA);
    cfg.trials = (ftdataICA.trialinfo==129);
    tl9 = ft_timelockanalysis(cfg, ftdataICA);
    dERP = tl7.avg-tl9.avg;
    dERN_66 = abs(dERP(lag_dERN)); % <--- THIS IS THE FINAL RESULT
end

% 100% Congruency
if sum(ftdataICA.trialinfo==137)==0 || sum(ftdataICA.trialinfo==139)==0
    dERN_100 = 0;
else
    cfg.trials = (ftdataICA.trialinfo==137);
    tl7 = ft_timelockanalysis(cfg, ftdataICA);
    cfg.trials = (ftdataICA.trialinfo==139);
    tl9 = ft_timelockanalysis(cfg, ftdataICA);
    dERP = tl7.avg-tl9.avg;
    dERN_100 = abs(dERP(lag_dERN)); % <--- THIS IS THE FINAL RESULT
end

%% Save results:

% dERN
dERN_output = [dERN dERN_0 dERN_33 dERN_66 dERN_100];
dlmwrite([cnfg.outpath 'dERN.txt'], dERN_output, 'delimiter', '\n');

% Summary
clear summary_txt
summary_txt{1} = ['Bad channels: ' BAD_CH{:}];
summary_txt{2} = ['Number of trials removed: ' num2str(BAD_TRIAL)];
summary_txt{3} = ['Selected Component: ' ftdataICA.label{ICoI}];
fid = fopen([cnfg.outpath 'summary.txt'], 'w');
for i = 1:length(summary_txt)
    fprintf(fid, '%s\n', summary_txt{i});
end
fclose(fid);





