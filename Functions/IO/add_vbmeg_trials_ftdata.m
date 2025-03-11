function ftdata_trials = add_vbmeg_trials_ftdata(cnfg,ftdata)

% Read from a VMBEG file some trials and create the ft structure
% accordingly to them
%
% cfg.filetrialname
% cfg.stimdef - in seconds [e.g. -0.2 0.5]
% cfg.trig - define the trigger to use, as there are several channels in
% bexp_ext
%
% Victor J. Lopez-Madrona
% 10/03/2025


%% Initialization
if ~isfield(cnfg,'trig'), cnfg.trig=1; end

warning('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
warning('Code in development. So far, only works with one trigger')
warning('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Pipeline
load(cnfg.filetrialname,'bexp_ext')
Ntrig = length(cnfg.trig);

for trigi = 1:Ntrig
    trig_aux = bexp_ext(cnfg.trig(trigi),:);
    % The th is the middle point between max and min
    th = max(trig_aux) - (max(trig_aux) - min(trig_aux)) / 2;
    trig_aux(bexp_ext(cnfg.trig(trigi),:)>=th) = 1;
    trig_aux(bexp_ext(cnfg.trig(trigi),:)<th) = 0;
    
    [~,p]=find(diff(trig_aux)==1);
    
    % Create trial configuration for fieldtrip
    Fs = ftdata.fsample; % Sampling frequency from FieldTrip data
    onset_samples = p; % Onset times in samples
    pre_samples = round(cnfg.stimdef(1)*Fs); % Time before onset in samples
    post_samples = round(cnfg.stimdef(2)*Fs); % Time after onset in samples

    % Define trials
    trl = zeros(length(onset_samples), 3);
    trl(:, 1) = onset_samples + pre_samples;  % Start sample
    trl(:, 2) = onset_samples + post_samples; % End sample
    trl(:, 3) = pre_samples;                  % Trial onset after sample starts 

    cfg = [];
    cfg.trl = trl;  % Assign the trial matrix
    ftdata_trials = ft_redefinetrial(cfg, ftdata);
    ftdata_trials.trialinfo = ones(length(trl),1)*trigi;
end
    
    
