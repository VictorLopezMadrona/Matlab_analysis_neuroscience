
function phase_profile = compute_phase_profile(x,freq,Fs)

%% Obtain the averaged shape of the cycles at one specific frequency
% For a given signal x and a frequency, the functions identifies all the
% cycles, resamples them to have the same duration and then averages
% them.
%
% This code is much faster than ITPC_MEG_SEEG 
%
% Syntax:
%    phase_profile = compute_phase_profile(x,freq,Fs)
%
% Inputs:
%   x    - signal (only one channel)
%   freq - [fmin fmax] range with the frequencies of interest. Ex: 7-9 Hz
%   Fs   - sampling frequency
%
% Outputs:
%   phase_profile - averaged cycle
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Jul. 2025; Last revision: 17-Jul-2025

%% Find the  cycles

% Filter the signal
xf = eegfilt(x,Fs,freq(1),freq(2));

%Find the phases with hilbert
ph = angle(hilbert(xf));

%Find starting point of each cycle
[~,start_cycle] = findpeaks(ph);

%Measure the duration (in time samples) of each cycle.
dur_cycle = diff(start_cycle);

%Select only those cycles with a duration in the range of the Frequency of
%interest (we accept 10% of variability)
min_dur = floor( Fs/freq(2) - 0.1*Fs/freq(2)); %minimal duration in samples
max_dur = floor( Fs/freq(1) + 0.1*Fs/freq(1)); %maximal duration in samples
select_cycle = find(dur_cycle>min_dur & dur_cycle<max_dur); %these are the cycles we select

% Determine the averaged duration of the cycles. It will be the cycle
% duration of the middle frequency. For example, if the freq range is 8-12 Hz,
% we select a common duration of 100ms (10 Hz) for all cycles
fc = freq(1)+diff(freq)/2;
mean_dur = round(Fs/fc);

% Loop across all identified cycles. We know the starting and endind point
% (start_cycle). So we are coming back to the original raw signal, select
% the time interval of each cycle, resample them to have the same
% duration (mean_dur) and average them.
phase_profile = zeros(1,mean_dur);
for i=select_cycle
    tini = start_cycle(i); %initial point of the cycle
    tend = start_cycle(i+1); %end point
    cycle_resampled = resample(x(tini+1:tend), mean_dur, length(tini+1:tend));
    phase_profile = phase_profile + cycle_resampled;
end
phase_profile = phase_profile/length(select_cycle);







