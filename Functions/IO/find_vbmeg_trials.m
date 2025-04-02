function trl = find_vbmeg_trials(cnfg,ftdata)

% Read from a VMBEG file some trials and return the time points.
% We do not create the trials here, just obtain the time information
%
% cfg.filetrialname
% cfg.stimdef - in seconds [e.g. -0.2 0.5]
% cfg.trig - define the trigger to use, as there are several channels in
% bexp_ext
%
% Victor J. Lopez-Madrona
% 02/04/2025


%% Initialization
if ~isfield(cnfg,'trig'), cnfg.trig=1; end

warning('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
warning('%% Code in development. So far, only works with one trigger %%')
warning('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Pipeline
load(cnfg.filetrialname,'bexp_ext','MEGinfo')
if ftdata.fsample ~= MEGinfo.SampleFreq
    error('The Sampling Frequency is different. Run this code before downsampling')
end

Ntrig = length(cnfg.trig);

for trigi = 1:Ntrig
    trig_aux = bexp_ext(cnfg.trig(trigi),:);
    % The th is the middle point between max and min
    th = max(trig_aux) - (max(trig_aux) - min(trig_aux)) / 2;
    trig_aux(bexp_ext(cnfg.trig(trigi),:)>=th) = 1;
    trig_aux(bexp_ext(cnfg.trig(trigi),:)<th) = 0;
    
    [~,p]=find(diff(trig_aux)==1);
    
    % Define trials
    Fs = ftdata.fsample; % Sampling frequency from FieldTrip data
    pre_samples = round(cnfg.stimdef(1)*Fs); % Samples before onset in samples
    post_samples = round(cnfg.stimdef(2)*Fs); % Samples after onset in samples
    
    trl = zeros(length(p), 3);
    trl(:, 1) = ftdata.time{1}(p+pre_samples);  % Start sample
    trl(:, 2) = ftdata.time{1}(p+post_samples); % End sample
    trl(:, 3) = ftdata.time{1}(p);  
end