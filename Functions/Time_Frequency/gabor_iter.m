function [gabor_mean,tt,ff] = gabor_iter(cnfg,ftdata)

%% From continuos data (FieldTrip or matrix), compute Gabor and TF response.
%
% The same as tf_gabor, but it is computed for each channel iteratively,
% which helps if not enough memory.
% Nothing is saved with this function.
%
% Syntax:
%    gabor_mean = gabor_iter(cnfg,ftdata)
%
% Inputs:
%    cfg - Structure of parameters:
%
%       stimdef    - [start end] Trial in seconds. Ex [-0.2 0.5]
%       time_trial - vector with each trial onset in seconds
%       M          - Time length for Gabor im samples (Def 128)
%       a          - Gabor resolution (Def M/16)
%       freqlim    - Frequencies of interest. Ex [5 100]
%       baseline   - [tmin tmax]; 'all' (def); 'none'
%
%   data - Data in format field trip or matrix.
%
% %%%%%%%%%%% CASE 'data' is a matrix %%%%%%%%%%%%%%%%
% This code has been prepared to work with FieldTrip structures. However,
% it can be used with any matrix of data by defining some extra parameters:
%
%   data        - [Nch x Nsamples]
%   cfg.time    - vector of Nsamples with the timestamp of each sample
%   cfg.Fs      - Sampling rate
%   cfg.label   - [Optional] cell of strings with the name of each channel
%
% Outputs:
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug. 2025; Last revision: 15-Apr-2025


%% Initialization

if ~isfield(cnfg,'M'), cnfg.M = 128; end
if ~isfield(cnfg,'a'), cnfg.a = cnfg.M/16; end
if ~isfield(cnfg,'baseline'), cnfg.baseline = 'all'; end
%if ~isfield(cnfg,'freqlim'), cnfg.freqlim = [0 ftdata.fsample/2]; end

%%% Case ftdata is a matrix
if ~isstruct(ftdata)
    [Nch,Nsamples] = size(ftdata);
    disp(['Input data is a matrix with ' num2str(Nch) ' channels and '...
        num2str(Nsamples) ' samples.'])
    if ~isfield(cnfg,'time') 
        error('Parameter ''time'' is mandatory when working with matrix')
    end
    if ~isfield(cnfg,'Fs') 
        error('Parameter ''Fs'' is mandatory when working with matrix')
    end
    if ~isfield(cnfg,'label'), cnfg.label = []; end
    
    data = ftdata;
    clear ftdata
    ftdata.trial{1} = data;
    ftdata.fsample  = cnfg.Fs;
    ftdata.time{1}  = cnfg.time;
    ftdata.label    = cnfg.label;
end

%% COMPUTE GABOR

% Parameters
M = cnfg.M; %Time length
a = cnfg.a; % I divide the Fs by this 'a' value
time_trial = cnfg.time_trial;

%%%% GABOR TRANSFORM %%%%
disp('Computing GABOR Transform...')

Fs = ftdata.fsample;
Nch = size(ftdata.trial{1},1);
chi = 0;
wb = waitbar(0,['Computing GABOR: channels ' num2str(chi) '/' num2str(Nch)]);

for chi=1:Nch
    
    y = ftdata.trial{1}(chi,:)';
    g = gabwin({'tight', 'hann'}, a, M);
    c = dgtreal(y, g, a, M);
    ff = linspace(0,Fs/2,size(c,1));
    
    %%% Correct power - NO BASELINE SELECTION FOR NOW
    % Case zscore each frequency
    if strcmp(cnfg.baseline,'all')
        pow_corrected = [];
        for fi=1:length(ff)
            pow_corrected(fi,:) = zscore(abs(c(fi,:)));
        end
        
        %No correction
    else
        pow_corrected = abs(c);
    end
    
    % Create trials from continuous GABOR. Here I need time_trial
    Ntrial = length(time_trial);
    Fs_gabor = round(ftdata.fsample/a);
    time_gabor = linspace(ftdata.time{1}(1),ftdata.time{1}(end),size(pow_corrected,2));
    % Find closest index and value for each y
    onset_samples = zeros(1,Ntrial);
    for tri = 1:Ntrial
        [~, onset_samples(tri)] = min(abs(time_gabor - time_trial(tri)));  % index of minimum difference
    end
    pre_samples   = round(cnfg.stimdef(1)*Fs_gabor); % Samples before onset in samples
    post_samples  = round(cnfg.stimdef(2)*Fs_gabor); % Samples after onset in samples
    % Define trials
    trl_ini = onset_samples + pre_samples;  % Start sample
    trl_end = onset_samples + post_samples; % End sample
    pow_trial_zscore = [];
    for tri=1:Ntrial
        pow_trial_zscore(:,:,tri) = pow_corrected(:,trl_ini(tri):trl_end(tri));
    end
    
    gabor_mean(:,:,chi) = mean(pow_trial_zscore,3);
    
    wb = waitbar(chi/Nch,wb,['Computing GABOR: channels ' num2str(chi) '/' num2str(Nch)]);
end
tt = linspace(cnfg.stimdef(1),cnfg.stimdef(2),size(gabor_mean,2));
ff = linspace(0,Fs/2,size(c,1));
close(wb)




