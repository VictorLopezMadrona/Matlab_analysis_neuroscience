
function ftdata_trial = create_time_trials(cnfg,ftdata)

% Create trial from a vector with the time points (time_trial)
%

time_trial = cnfg.time_trial;

Ntrial = length(time_trial);
Fs = ftdata.fsample;
% Find closest index and value for each y
onset_samples = zeros(1,Ntrial);
for tri = 1:Ntrial
    [~, onset_samples(tri)] = min(abs(ftdata.time{1} - time_trial(tri)));  % index of minimum difference
end
pre_samples   = round(cnfg.stimdef(1)*Fs); % Samples before onset in samples
post_samples  = round(cnfg.stimdef(2)*Fs); % Samples after onset in samples

% Define trials
trl = zeros(length(onset_samples), 3);
trl(:, 1) = onset_samples + pre_samples;  % Start sample
trl(:, 2) = onset_samples + post_samples; % End sample
trl(:, 3) = pre_samples;                  % Trial onset after sample starts

cfg = [];
cfg.trl = trl;  % Assign the trial matrix
ftdata_trial = ft_redefinetrial(cfg, ftdata);
ftdata_trial.trialinfo = ones(length(trl),1);


