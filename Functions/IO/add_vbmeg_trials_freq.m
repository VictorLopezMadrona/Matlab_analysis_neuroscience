function freq_trials = add_vbmeg_trials_freq(cnfg,freq)

% Read from a VMBEG file some trials and divide the time-freq into trials
%
% cfg.filetrialname
% cfg.stimdef - in seconds [e.g. -0.2 0.5]
% cfg.trig - define the trigger to use, as there are several channels in
% cfg.Fs
% bexp_ext
%
% Victor J. Lopez-Madrona
% 13/03/2025


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
    Fs = cnfg.Fs; % Sampling frequency from FieldTrip data
    onset_time = p/Fs;
 
    % Define trials
    trl = zeros(length(onset_time), 3);
    trl(:, 1) = onset_time + cnfg.stimdef(1);  % Start sample
    trl(:, 2) = onset_time + cnfg.stimdef(2); % End sample
end

% Initialize a new freq structure for trials
freq_trials = freq;
freq_trials.fourierspctrm = [];
freq_trials.powspctrm = [];
freq_trials.time = [];

Ntrial = size(trl,1);

for t = 1:Ntrial
    % Find the time indices corresponding to each trial
    t_start = freq.time >= trl(t,1);
    t_end = freq.time <= trl(t,2);
    
    % Extract the time range for this trial
    trial_idx = t_start & t_end;
    
    % Store the power spectrum for this trial
    if isfield(freq,'fourierspctrm')
        freq_trials.fourierspctrm(t,:,:,:) = freq.fourierspctrm(1,:,:,trial_idx);
    end
    if isfield(freq,'powspctrm')
        freq_trials.powspctrm(t,:,:,:) = freq.powspctrm(:,:,trial_idx);
    end
    
    % Store the corresponding time vector
    freq_trials.time{t} = freq.time(trial_idx);
end

% Adjust the dimension labels
freq_trials.dimord = 'rpt_chan_freq_time';  % Trials, Channels, Frequencies, Time

    
    








