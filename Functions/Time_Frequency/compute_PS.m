
function [PS,f,fbands] = compute_PS(x,nfft,overlap,Fs,bands)
% Function to cumpute the Power Spectrum using Fast Fourier Transform and a
% sliding window.
%
% USE: [PS,f,fbands] = compute_PS(x,nfft,overlap,Fs, band);
%
% INPUT:
%   x [Nch x samples]: input signal.
%   nfft: Number of frequencial points. Optional (default: length(x)).
%   overlap: Number of point that overlap each window. Optional (default: 0).
%   Fs: Sampling Frequency. Optional (default: 250 Hz).
%   band [Nbands x 2]: Freq bands to compute the averaged PS (default: 0).
%                      Example: band = [1 4; 4 8; 8 12];
%
% OUTPUT:
%   PS [Nch x nfft]: Module of the Power Spectrum for each channel.
%   f: vector with frequency.
%
% Víctor José López Madrona - 14/6/2018


if nargin == 1
    nfft = length(x);
    overlap = 0;
    Fs = 250;
    bands = 0;
elseif nargin == 2
    overlap = 0;
    Fs = 250;
    bands = 0;
elseif nargin == 3
    Fs = 250;
    bands = 0;
elseif nargin == 4
    bands = 0;
end

if overlap >= nfft || overlap < 0
    error('The overlaping parameter must be between 0 and nfft')
end
overlap = nfft-overlap; %Instead of overlaping, we use the displacement of the window

x = norm_signal(x); %In case of normalizate the signals...
[nch,samples] = size(x);
Nw = floor( (samples - nfft) / overlap) + 1;
PS = zeros(nch,nfft);
for n=1:nch
    for w=1:Nw
        PS(n,:) = PS(n,:) + abs(fft(x(n,(w-1)*overlap+1:(w-1)*overlap+nfft)));
    end
end
PS = PS/Nw;
PS = PS/nfft;
f = linspace(0,Fs,nfft);


%Compute average fft for each band
if bands == 0
    fbands=NaN;
else
    [Nbands,~] = size(bands);
    fbands = zeros(nch,Nbands);
    size(PS)
    for b=1:Nbands 
        fbands(:,b) = mean(PS(:,f>=bands(b,1) & f<=bands(b,2)),2);
    end
end

%% PLOT RESULTS
% Example code to plot the results of one of the channels.

% ch = 1; %Channel to plot
% figure, 
% subplot(1,2,1)
% plot(f,PS(ch,:));
% xlabel('Frequency (Hz)')
% ylabel('Power Spectrum (mV^2/Hz')
% subplot(1,2,2)
% bar(fbands(ch,:))
% ylabel('Power Spectrum (mV^2/Hz)')
% for b=1:Nbands
%     xtl{b}=[num2str(bands(b,1)) '-' num2str(bands(b,2))];
% end
% set(gca,'XTickLabels',xtl)

function x_out = norm_signal(x_in)

% Function to normalize a signal, making its mean zero and variance unit.
%
% USE:
%   x_out = norm_signal(x_in);
%
% INPUT:
%   x_in [NCh x Samples]: Matrix with our signals.
%
% OUTPUT:
%   x_out [NCh x Samples]: Matrix with the normalized signals.
%
% Víctor José López Madrona - 06/03/2017

[nch,~] = size(x_in);
x_out=x_in*0;
for n=1:nch
    x_out(n,:) = x_in(n,:)/std(x_in(n,:));
    x_out(n,:) = x_out(n,:) - mean(x_out(n,:));
end

%Normalizate using the mean and variance of the whole dataset.
clear x_out
x_out = x_in/std(x_in(:));
x_out = x_out - mean(x_out(:));



        

        



