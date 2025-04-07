
%% Example NMF on a simple simulation

%% Generate a toy model with beta increase or decrease

%%% Parameters %%%
fs = 512;              % Sampling frequency (Hz)
duration = 1;          % Duration (s)
freq_burst = [15 25];  % Frequency band of the activity to include
A_burst = 3;           % Amplitude of the burst with respect to the noise
dur_burst = [0.3 0.5]; % Time window where the burst appears. The first point will be considered as trial onset
Ntrial = 100;          % Number of trials
% vector with the number of 'sensors' and the relative amplitude of the
% activity burst on each of them. 
sensor = [0.33 0.66 1 0.66 0.33];
beta_decrease = false; %Turn it 'true' for a beta supression

N = fs * duration;    % Number of samples
t = (0:N-1)/fs;       % Time vector
signal_trial = zeros(length(sensor),N,Ntrial);
signal_cont = [];
signal_beta = [];
w = (hamming(N).^0)'; % Generate Hamming window for each trial
w = repmat(w,length(sensor),1);

%%% Noise %%%
% Generate 1/f noise in frequency domain for all the time. Then we can cut
% into trials. Better than noise for each trial
Ntot = N*Ntrial;
f = (0:Ntot-1)*(fs/Ntot);       % Frequency vector
amplitude = 1 ./ (f + 1); % 1/f scaling
phase = 2*pi*rand(1, Ntot);  % Random phase
spectrum = amplitude .* exp(1i * phase);
spectrum(Ntot/2+2:end) = conj(spectrum(Ntot/2:-1:2));  % Enforce conjugate symmetry

% Inverse FFT to obtain time-domain signal
pink_noise_tot = real(ifft(spectrum));
pink_noise_tot = zscore(pink_noise_tot); %Normalize
pink_noise_ind = zeros(length(sensor),length(pink_noise_tot));

%%% Create time-course signal
for tri=1:Ntrial

    %Select one segment of the pink_noise
    pink_noise = pink_noise_tot((tri-1)*N+1:tri*N);
    
    %%% Activity burst %%%
    % Generate beta activity by filtering noise 
    white_noise = randn(1,length(pink_noise));
    activity_burst = eegfilt(white_noise,fs,freq_burst(1),freq_burst(2));

    % Do a Hanning window on the time of interest to apply to my activity burst
    start_time = dur_burst(1);
    end_time = dur_burst(2);
    % Indices for the window
    win_idx = t >= start_time & t <= end_time;
    win_len = sum(win_idx);  % Number of points in the window
    % Create Hanning window and place it in full-length array
    hann_window = zeros(size(t));
    hann_window(win_idx) = hann(win_len)';  % Transpose to match shape
    if beta_decrease
        hann_window = 1-hann_window;
    end

    %%% CREASE SIGNAL IN TRIALS %%%
    % Pink noise + (Amplitude*activity_burst).*Hanning window
    for si = 1:length(sensor)   
        % Each sensor will have the same pink noise + an independent pink noise
        f = (0:N-1)*(fs/N);       % Frequency vector
        amplitude = 1 ./ (f + 1); % 1/f scaling
        phase = 2*pi*rand(1, N);  % Random phase
        spectrum = amplitude .* exp(1i * phase);
        spectrum(N/2+2:end) = conj(spectrum(N/2:-1:2));
        pink_noise_sensor = real(ifft(spectrum));
        pink_noise_sensor = zscore(pink_noise_sensor);
        pink_noise_ind(si,(tri-1)*N+1:tri*N) = pink_noise_sensor;
        
        signal_trial(si,:,tri) = pink_noise + 0.35*pink_noise_sensor...
            + (sensor(si)*A_burst*activity_burst).*hann_window;
        
    end
        %signal_trial(tri,:) = pink_noise + (A_burst*activity_burst).*hann_window;
    
    %%% SIGNAL CONT %%%
    signal_cont = [signal_cont signal_trial(:,:,tri).*w];
    signal_beta = [signal_beta activity_burst.*hann_window];
    
    % % Plot time-domain signal
    % figure;
    % subplot(2,1,1);
    % plot(t, signal);
    % xlabel('Time (s)');
    % ylabel('Amplitude');
    % title('Time-Domain Signal: 1/f Noise with Beta Burst');
    % 
    % % Plot Power Spectral Density using Welch’s method
    % subplot(2,1,2);
    % pwelch(signal, [], [], [], fs);
    % title('Power Spectral Density (Welch Method)');
end

%%% All the data
signal_cont; %Continuous data
pink_noise_tot; %General noise to all channels
signal_beta; %beta activity for each channel. Related to 'sensor', which is the amplitude at each sensor
pink_noise_ind; %Individual noise for each channel

%% Compute NMF

% Obtain a vector with the starting point of each trial in seconds
% The starting point is the first value of 'dur_burst'
time_trial = dur_burst(1) : duration : (duration*(Ntrial-1)+dur_burst(1));
% The trial limits are defined by trial duration. 
% From 0 to dur_burst(1) is pre-trial. From dur_burst(1) to duration is
% post-trial
stimdef = [-dur_burst(1) duration-dur_burst(1)];

fs = 512;              % Sampling frequency (Hz)
duration = 1;          % Duration (s)
freq_burst = [15 25];  % Frequency band of the activity to include
A_burst = 3;           % Amplitude of the burst with respect to the noise
dur_burst = [0.3 0.5]; % Time window where the burst appears. The first point will be considered as trial onset
Ntrial = 100;  

%%% Set the parameters and run the function
data = [];
data = signal_cont; %This is my input data

cfg = [];
cfg.time       = t; % Time vector
cfg.Fs         = fs;
cfg.time_trial = time_trial;
cfg.stimdef    = stimdef; %Time onset of each trial
cfg.doplot     = true; %Plot the result
cfg.dosave     = false; %Do not save the result
cfg.freqlim    = [5 100]; %Frequencies where the NMF is computed
cfg.M          = 128; %Length of time window for time-frequency
cfg.a          = cfg.M/8; %Step of the time window
cfg.stats      = 'std'; %Select outlayers based on 2x standard deviation
cfg.alpha      = 0.01;

[Wtf,H,W] = tf_nmf(cfg,data);










