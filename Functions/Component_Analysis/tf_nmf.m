function [Wtf,H,W] = tf_nmf(cnfg,ftdata)

%% From continuos data in FieldTrip format, compute Gabor, trials and NMF.
% It will do the trial average before the NMF.
% Time and Frequency are vectorize into a single dimension
% The second dimension is space (sensor)
% Therefore, it will compute:
%
%   X ~= WH
%
% - where X is the Gabor transform with time-freq in one dimension and
% sensors (space) in the second dimension
% - W are the tf distribution of each component
% - H are the spatial distribution of each component
%
% Syntax:
%    [Wtf,H,W] = tf_nmf(cfg,ftdata)
%
% Inputs:
%    cfg - Structure of parameters:
%
%       Ncomp      - Number of components 'k'. It can be a value or it can
%                    be estimated ('all','velicer')
%       stimdef    - [start end] Trial in seconds. Ex [-0.2 0.5]
%       time_trial - vector with each trial onset in seconds
%       M          - Time length for Gabor im samples (Def 128)
%       a          - Gabor resolution (Def M/16)
%
%       dosave   - logical. True/false save/not save the results
%       outpath  - string. Path to save the results if dosave=true
%       plotfig  - logical. Plot the figure with the results.
%       infosave - string to include in the saved filed
%
%   ftdata - Data in format field trip
%
% Outputs:
%    Wtf - W matrix reshaped to be 2D (time x freq)
%    H   - Spatial distribution
%    W   - TimeFreq distribution (vectorized)
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug. 2025; Last revision: 02-Apr-2025


%% Initialization

if ~isfield(cnfg,'Ncomp'), cnfg.Ncomp = 'velicer'; end
if ~isfield(cnfg,'M'), cnfg.M = 128; end
if ~isfield(cnfg,'a'), cnfg.a = cnfg.M/16; end
    
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'plotfig'), cnfg.plotfig=true; end
if ~isfield(cnfg,'outpath') && cnfg.dosave
    error('Outpath has not been specified to save the results'), end

%% COMPUTE GABOR

% Parameters
M = cnfg.M; %Time length
a = cnfg.a; % I divide the Fs by this 'a' value
time_trial = cnfg.time_trial;

%%%% GABOR TRANSFORM %%%%
disp('Computing GABOR Transform...')
y = ftdata.trial{1}';
Fs = ftdata.fsample;
g = gabwin({'tight', 'hann'}, a, M);
c = dgtreal(y, g, a, M);
ff = linspace(0,Fs/2,size(c,1));
pow_corrected = [];
for si=1:size(c,3)
    for fi=1:length(ff)
        pow_corrected(fi,:,si) = zscore(abs(c(fi,:,si)));
    end
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
    pow_trial_zscore(:,:,:,tri) = pow_corrected(:,trl_ini(tri):trl_end(tri),:);
end

pow_mean = mean(pow_trial_zscore,4);

%% NMF

%%% Vectorize
[T, F, S] = size(pow_mean);
pow_vec = reshape(pow_mean, [T*F, S]); 
V = pow_vec - min(pow_vec(:));

%%% Number of components
if isnumeric(cnfg.Ncomp)
    k = cnfg.Ncomp;
elseif strcmp(cnfg.Ncomp,'all')
    k = size(V,2);
elseif strcmp(cnfg.Ncomp,'velicer')
    [~, k] = velicer_map_oconnor(V);
end

% SVD-based initialization for stability
[U,D,S] = svd(V,"econ");
S = S';
U = U*D;
W0 = abs(U(:,1:k));
H0 = abs(S(1:k,:));

[W, H] = nnmf(V, k, 'w0', W0, 'h0', H0, 'algorithm', 'mult');

%%% Reshape to original dimensions
[T, F, S] = size(pow_mean);
Wtf = reshape(W, [T, F, k]);   % Back to time × freq × sensors







