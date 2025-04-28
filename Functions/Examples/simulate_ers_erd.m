
function [data,ERS,tt] = simulate_ers_erd(cnfg)

%% Simulate realistic brain data with trial-based induced increases of 
% activity at a given frequency band.
% The ERD is generated as a general increase all the time instead of in the
% interval of interest
%
% By default, it uses pink noise as background noise. However, it accepts a
% signal as input. In that case, it generates a background noise with the
% same probability distribution and power spectrum (or very similar). Note
% that the output will be different each time the code is executed, even if
% the input data is the same.
%
% USE:
%   data = simulate_ers_erd(cfg);
%
% INPUT:
%   cfg.data - [Optional] It is possible to introduce a previous
%              time-course. If so, the simulation will have the same Power 
%              Spectrum as this data.
%
%      .Fs   - Sampling frequency of the input data OR for the pink noise (Def: 500)
%
%      .Ntrial    - Number of trials (Def = 100)
%      .trial_dur - Duration of each trial in seconds (Def = 2)
%      .time_band - Time interval for the induced ERS 
%                   (Def: [0.5 1]) The ERS appear between 0.5 and 1 s
%                   after trial onset
%      .freq_band - Frequency band for the induced ERS in Hz (Def: [70 - 90])
%      .amplitude - Amplitude of the ERS (Def: 1)
%      .plot_tf   - Compute and plot the time-freq. (Def = true)
%      .mode      - 'ers' or 'erd' (Def: 'ers')
%
% OUTPUT:
%   data - Simulated time-course with all the trials. 
%          data(1) corresponds to the first trial onset. 
%   ERS  - ERS (or ERD) time course
%   time-data - vector with the time points associated to 'data'
%
%   It also plots the time-frequency of the response using Gabor
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Apr. 2025; Last revision: 25-Apr-2025

%% Parameters

if nargin == 0, cnfg = []; end
if isfield(cnfg,'data') && ~isfield(cnfg,'Fs')
    error('You need to define the sampling frequency of the input data'), end
if ~isfield(cnfg,'Ntrial'), cnfg.Ntrial = 100; end
if ~isfield(cnfg,'trial_dur'), cnfg.trial_dur = 2; end
if ~isfield(cnfg,'time_band'), cnfg.time_band = [0.5 1]; end
if ~isfield(cnfg,'freq_band'), cnfg.freq_band = [70 90]; end
if ~isfield(cnfg,'amplitude'), cnfg.amplitude = 1; end
if ~isfield(cnfg,'Fs'), cnfg.Fs = 500; end
if ~isfield(cnfg,'plot_tf'), cnfg.plot_tf = true; end
if ~isfield(cnfg,'mode'), cnfg.mode = 'ers'; end

Fs         = cnfg.Fs;        % Sampling frequency (Hz)
duration   = cnfg.trial_dur; % Duration (s)
freq_burst = cnfg.freq_band; % Frequency band of the activity to include
dur_burst  = cnfg.time_band; % Time window where the burst appears. The first point will be considered as trial onset
Ntrial     = cnfg.Ntrial;    % Number of trials
A_burst    = cnfg.amplitude; % Amplitude of the burst with respect to the noise

%% Generate background noise

if ~isfield(cnfg,'data')
    disp('Computing simulation from pink noise')
    
    N = Fs * duration; % Number of samples per trial

    % Generate 1/f noise in frequency domain for all the time.
    Ntot = N*Ntrial;
    f = (0:Ntot-1)*(Fs/Ntot);       % Frequency vector
    amplitude = 1 ./ (f + 1); % 1/f scaling
    phase = 2*pi*rand(1, Ntot);  % Random phase
    spectrum = amplitude .* exp(1i * phase);
    spectrum(Ntot/2+2:end) = conj(spectrum(Ntot/2:-1:2));  % Enforce conjugate symmetry

    % Inverse FFT to obtain time-domain signal
    background_noise = real(ifft(spectrum));
    background_noise = zscore(background_noise); %Normalize

else
    disp('Computing simulation from input data')
    disp(['Sampling frequency of input data is: ' num2str(Fs)])
    
    N = Fs * duration; % Number of samples per trial
    Ntot = N*Ntrial; % Total number of samples
    
    %If the input data is too short, we mirror it until we have enough length
    x = cnfg.data;
    Nw = ceil(Ntot,length(x)); 
    x_aux = x;
    for wi=1:Nw
        x = x(end:-1:end); %mirror the signal
        x_aux = [x_aux x]; %concatenate the mirror signal
    end
    
    % Generate surrogated data with the same prob dist and power spectrum
    background_noise = generate_iAAFT(x_aux);
    background_noise = background_noise(1:Ntot);
end

%% Generate and include ERS
    
% Do a Hanning window on the time of interest to apply to my activity burst
tt = (0:N-1)/Fs; % Time vector
start_time = dur_burst(1);
end_time = dur_burst(2);
% Indices for the window
win_idx = tt >= start_time & tt <= end_time;
win_len = sum(win_idx);  % Number of points in the window
% Create Hanning window and place it in full-length array
hann_trial = zeros(size(tt));
hw = hann(win_len).^0.2;
hann_trial(win_idx) = hw';  % Transpose to match shape

if strcmp(cnfg.mode,'ers')
    hann_all = repmat(hann_trial,1,Ntrial);
elseif strcmp(cnfg.mode,'erd')
    hann_all = repmat(hann_trial,1,Ntrial);
    hann_all = 1 - hann_all;
else
    error('cfg.mode was not recognized')
end
    

% Generate oscillation
cfg = [];
cfg.Fs = cnfg.Fs;      
cfg.dur = Ntrial*duration+1;      
cfg.freq_band = cnfg.freq_band; 
ERS = generate_oscillation(cfg);
% Filter the signal to avoid sharp transients
[b, a] = butter(2, [freq_burst(1) freq_burst(2)]/(Fs/2), 'bandpass');  % 2nd-order Butterworth
ERS = filtfilt(b, a, ERS);       % Zero-phase filtering same noise
% Normalize the oscillation to have the same mean and std than the noise
mu = mean(background_noise);
sigma = std(background_noise);
ERS = (ERS/std(ERS))*sigma;
ERS = ERS-mean(ERS)+mu;
ERS = ERS(1:length(hann_all));

% Filter the noise around the frequency band of interest
%[b, a] = butter(2, [freq_burst(1) freq_burst(2)]/(Fs/2), 'bandpass');  % 2nd-order Butterworth
%ERS = filtfilt(b, a, background_noise);       % Zero-phase filtering same noise
%ERS = filtfilt(b, a, background_noise);       % Reversed noise
%ERS = filtfilt(b, a, randn(1,Ntot));       % White noise   
%ERS = eegfilt(randn(1,Ntot), Fs, freq_burst(1),0); % 
%ERS = eegfilt(ERS, Fs, 0,freq_burst(2)); % 

% The data is the backgound noise + the filter data after the hanning
% window and ultiply by the desired amplitude
ERS = A_burst.*ERS.*hann_all;
data = background_noise + ERS;
time_data = linspace(0,Ntrial*duration,length(data));

%% Do time-frequency

if cnfg.plot_tf

    y = data';
    
    %%%% GABOR TRANSFORM %%%%
    T = length(y(:,1));
    M = 128; %Time length
    a = M/32;
    g = gabwin({'tight', 'hann'}, a, M);
    c = dgtreal(y, g, a, M);
    
    tt = linspace(0,time_data(end),size(c,2));
    ff = linspace(0,Fs/2,size(c,1));
    %figure, imagesc(tt,ff,abs(c))

    %%% z-score in all the time-course %%%
    pow_corrected = [];
    for fi=1:length(ff)
        pow_corrected(fi,:) = zscore(abs(c(fi,:)));
    end
    figure, imagesc(tt,ff,pow_corrected)
    axis('xy')
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('Continues TF - zscore each freq')
 
    % Create trial corrected
    pow_trial = [];
    tr_dur = floor(size(c,2)/Ntrial);
    for tri=1:Ntrial
        time_ini = (tri-1)*duration;
        tr_ini = find(abs(tt-time_ini)>0,1);
        pow_trial(:,:,tri) = pow_corrected(:,tr_ini:tr_ini+tr_dur-1);
        %pow_trial_zscore(:,:,tri) = pow_corrected(:,tr_dur*(tri-1)+1:tr_dur*tri);
    end
    figure,
    tt = linspace(0,duration,size(pow_trial,2));
    imagesc(tt,ff,squeeze(mean(pow_trial,3)))
    axis('xy')
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('Averaged TF - zscore each freq')
    
    %%% Baseline correction %%%
    gab_trial = [];
    tr_dur = floor(size(c,2)/Ntrial);
    for tri=1:Ntrial
        time_ini = (tri-1)*duration;
        tr_ini = find(abs(tt-time_ini)>0,1);
        gab_trial(:,:,tri) = c(:,tr_ini:tr_ini+tr_dur-1);
    end
    
    baseline_time = 1:5; % First 5 samples
    pow_trial = [];
    for fi=1:length(ff)
        for tri=1:Ntrial
            baseline = mean(abs(gab_trial(fi,baseline_time,tri)));
            pow_trial(fi,:,tri) = (abs(gab_trial(fi,:,tri)) - baseline)/baseline;
        end
    end
    figure, 
    tt = linspace(0,duration,size(pow_trial,2));
    imagesc(tt,ff,squeeze(mean(pow_trial,3)))
    axis('xy')
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title('Averaged TF - baseline corrected')
end



