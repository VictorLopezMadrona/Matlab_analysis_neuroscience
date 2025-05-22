
function [signal, t] = generate_oscillation(cnfg)

%% Generates an oscillatory signal with a given bandwidth.
% The duration of each cycle is obtained from a gaussian distribution
% centered at the central frequency and with a given bandwidth as standard
% deviation
%
% USE:
%   x = generate_oscillation(cfg);
%
% INPUT:
%   cfg.Fs        - Sampling frequency
%      .dur       - Total duration in seconds.
%      .freq_band - Frequency band of the oscillation
%
% OUTPUT:
%   x - 
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Apr. 2025; Last revision: 28-Apr-2025

Fs = cnfg.Fs; % Sampling frequency (Hz)
T_total = cnfg.dur; % Total duration (seconds)
t = 0:1/Fs:T_total; % Time vector

signal = zeros(size(t)); % Pre-allocate signal
current_time = 0;

%dur_cycle = 1./cnfg.freq_band;
%diff_dur = dur_cycle(1)-dur_cycle(2);

BW = cnfg.freq_band(2)-cnfg.freq_band(1);
Fc = cnfg.freq_band(1) + BW/2;
dur_cycle = 1./cnfg.freq_band;

%mu = 1/Fc;    % Mean cycle duration in seconds
%sigma = (1/cnfg.freq_band(1))/4; % Standard deviation in seconds

mu = Fc;
sigma = BW/4;

while current_time < T_total
    % Random cycle duration from normal distribution
    freq_cycle = normrnd(mu, sigma);
    cycle_duration = 1/freq_cycle;
    %cycle_duration = normrnd(mu, sigma);
    
    % Keep cycle duration within the BW
    cycle_duration = max(min(cycle_duration, dur_cycle(1)), dur_cycle(2));
    
    % Corresponding frequency
    freq = 1 / cycle_duration;
    
    % Time vector for this cycle
    t_cycle = current_time:1/Fs:(current_time+cycle_duration);
    if t_cycle(end) > T_total
        t_cycle = current_time:1/Fs:T_total;
    end
    
    % Generate one cycle
    signal_cycle = sin(2*pi*freq*(t_cycle - current_time));
    
    % Insert cycle into full signal
    idx_start = round(current_time*Fs) + 1;
    idx_end = idx_start + length(signal_cycle) - 1;
    signal(idx_start:idx_end) = signal_cycle;
    
    % Update current time
    current_time = current_time + cycle_duration;
end

% Plot
% figure,
% plot(t, signal)
% xlabel('Time (s)')
% ylabel('Amplitude')
% title('Oscillation with Gaussian Cycle Durations')


