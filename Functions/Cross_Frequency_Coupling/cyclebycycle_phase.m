
function [x_phase, phase_symmetry] = cyclebycycle_phase(x,f_phase,Fs,BW)

%% Measure the phase of one signal based on the "cycle-by-cycle"
% analysis by Cole and Voytek (Cycle-by-cycle analysis of neural oscillations, 2018) 

% It computes also "theta symmetry" following (Hippocampal theta bursting 
% and waveform shape reflect CA1 spiking patterns, 2018)
%
% USE:
%   [x_phase, phase_symmetry] = cyclebycycle_phase(x,f_phase,Fs,BW)
%
% INPUT:
%
% OUTPUT:
%
% See also: computeCFC

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Apr. 2018; Last revision: 21-Apr-2025

if nargin==1
    f_phase = 8;
    Fs = 625;
    BW = 1;
elseif nargin==2
    Fs = 625;
    BW = 1;
elseif nargin==3
    BW = 1;
end

xB = eegfilt(x,Fs,0,25);
xB = eegfilt(xB,Fs,1,0);
xC = eegfilt(x,Fs,0,f_phase+BW);
xC = eegfilt(xC,Fs,f_phase-BW,0);
%xhil = angle(hilbert(xC));
[~,x0up]=find(diff(xC>0)==1); x0up=x0up+1;
[~,x0dw]=find(diff(xC>0)==-1); x0dw=x0dw+1;
%[~,x0up]=find(diff(xhil>=-pi/2)==1); x0up=x0up+1;
%[~,x0dw]=find(diff(xhil>=pi/2)==1); x0dw=x0dw+1;

[~,d]=find(x0dw>x0up(1),1); %Ver si tengo que eliminar el último ciclo
if d>1, x0dw(1:d-1)=[]; end

% STEP D+E:
for l=1:min([length(x0dw) length(x0up)])-1
    [Pm(l),Ps(l)] = max(xB(x0up(l):x0dw(l))); Ps(l)=Ps(l)+x0up(l); %Maximum
    [Tm(l),Ts(l)] = min(xB(x0dw(l):x0up(l+1))); Ts(l)=Ts(l)+x0dw(l); %Minimum
    
    %Instead of finding when xB is halfway between peaks, we find the closest point 
    %in xB to x0dw (zero-crossing of step C) and y0dw (voltage halfway).
    y0dw(l)=Pm(l)-(Pm(l)-Tm(l))/2; %Middle point between max and min
    d=sqrt(((Ps(l):Ts(l))-x0dw(l)).^2 + ((xB(Ps(l):Ts(l)))-y0dw(l)).^2);
    [~,p]=min(d); C0dw(l)=p+Ps(l)-1;
    
    if l>1
    y0up(l)=Pm(l)-(Pm(l)-Tm(l-1))/2;
    d=sqrt(((Ts(l-1):Ps(l))-x0up(l)).^2 + ((xB(Ts(l-1):Ps(l)))-y0up(l)).^2);
    [~,p]=min(d); C0up(l-1)=p+Ts(l-1)-1;
    end
end
%C0dw=x0dw(1:end);
%C0up=x0up(2:end);

% STEP E:
% for l=1:length(Ps)-1
%     [~,C0dw(l)] = find(xB(Ps(l):Ts(l)) < Pm(l)-(Pm(l)-Tm(l))/2,1); C0dw(l)=C0dw(l)+Ps(l);
%     C0dwv(l)=Pm(l)-(Pm(l)-Tm(l))/2;
%     [~,C0up(l)] = find(xB(Ts(l):Ps(l+1)) > Pm(l+1)-(Pm(l+1)-Tm(l))/2,1); C0up(l)=C0up(l)+Ts(l);
%     C0upv(l)=Pm(l+1)-(Pm(l+1)-Tm(l))/2;
% end


x_phase=x*0;
for l=1:length(Ps)-1
    x_phase(Ps(l):C0dw(l)) = linspace(0,pi/2,length(Ps(l):C0dw(l)));
    x_phase(C0dw(l):Ts(l)-1) = linspace(pi/2,pi,length(C0dw(l):Ts(l)-1));
    x_phase(Ts(l):C0up(l)) = linspace(-pi,-pi/2,length(Ts(l):C0up(l)));
    x_phase(C0up(l):Ps(l+1)) = linspace(-pi/2,0,length(C0up(l):Ps(l+1)));
end

%symmetry
for l=2:length(Ps)-1
    rise_decay(l) = (Ps(l+1)-Ts(l)) / (Ps(l+1)-Ps(l));
    peak_trough(l) = (C0dw(l)-C0up(l-1)) / (C0up(l)-C0up(l-1));
end
    
phase_symmetry.rise_decay = rise_decay;
phase_symmetry.peak_trough = peak_trough;
phase_symmetry.peaks =  Ps;
phase_symmetry.trough =  Ts;
phase_symmetry.C0_down =  C0dw;
phase_symmetry.C0_up =  C0up;