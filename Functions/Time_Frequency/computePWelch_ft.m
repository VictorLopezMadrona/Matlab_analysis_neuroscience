function [Power,F]=computePWelch_ft(ftdata,cnfg)

%% Computes the Power Spectrum using PWelch
% It computes the PS for each trial and averages all of them
%
% USE: PS=computePS_ft(ftdata,cfg)

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug 21; Last revision: 24-Aug-2021

if nargin == 1
    cnfg=[];
end
if ~isfield(cnfg,'nfft'), cnfg.nfft=round(ftdata.fsample*4); end
if ~isfield(cnfg,'window'), cnfg.window=round(cnfg.nfft/2); end
if ~isfield(cnfg,'noverlap'), cnfg.noverlap=round(cnfg.window/2); end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if isfield(cnfg,'outpath')
    if ~strcmp(cnfg.outpath(end),'\')
        cnfg.outpath = [cnfg.outpath '\'];
    end
end
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end

%% PIPELINE


Fs=ftdata.fsample;
nfft=cnfg.nfft;
window=cnfg.window;
noverlap=cnfg.noverlap;
data = ftdata.trial{1};

[Power,F]   = pwelch(data',window,noverlap,nfft,Fs);

%% PLOT RESULTS

%%% PLOT POWER SPECTRUM ANALYSIS
nrows       = ceil(sqrt(length(ftdata.label)));
ncols       = ceil(length(ftdata.label)/nrows);
h=figure;
for ch=1:length(ftdata.label)
    G(ch) = subplot(nrows,ncols,ch);
    hold on,
    P = plot(F,Power(:,ch),'Linewidth',1.5);
    %P = plot(F,Power(:,ch).*F,'Linewidth',1.5); %Correct 1/f
    % Define action on subplot click and disble content selection
    set(P,'UserData',ch);
    set(P, 'HitTest', 'off'); % Disable content selection
    G(ch).ButtonDownFcn = @newFigure1; % Set click action
    
    xlim([1 49])
    %ylim([min_val max_val])
    title(ftdata.label{ch})
    
    set(G(ch),'xscale','log')
    xticks([1 2 4 8 16 32 64 128])

end

if cnfg.dosave
    if ~exist([cnfg.outpath '\Pwelch'],'dir'), mkdir([cnfg.outpath '\Pwelch']); end
    savefig(h,[cnfg.outpath '\Pwelch\Pwelch'])
    saveas(h,[cnfg.outpath '\Pwelch\Pwelch.png'])
    save([cnfg.outpath '\Pwelch\Pwelch_val'],'F','Power')
end

% Repeat the figure in log
for ch=1:length(ftdata.label)
    G(ch) = subplot(nrows,ncols,ch);
    xlim([1 150])
    set(G(ch),'xscale','log')
    set(G(ch),'yscale','log')
    xticks([1 2 4 8 16 32 64 128])
    xticklabels({'1' '2' '4' '8' '16' '32' '64' '128'})
end
if cnfg.dosave
    saveas(h,[cnfg.outpath '\Pwelch\Pwelch_log.png'])
end

% Aux figure for other parameters
for ch=1:length(ftdata.label)
    G(ch) = subplot(nrows,ncols,ch);
    xlim([2 40])
    set(G(ch),'xscale','log')
    set(G(ch),'yscale','linear')
    xticks([1 2 4 8 16 32 64 128])
    xticklabels({'1' '2' '4' '8' '16' '32' '64' '128'})
end

if ~cnfg.doplot
    close(h);
end
end

function newFigure1(h1,~)
%% Function to act on subplot click
%         Mouse click: Plots the selected subplot to a new figure
%  Ctrl + Mouse click: Delete subplot

switch get(gcf,'SelectionType')
    case 'normal'
        F = figure();
        copyobj(h1.Children,gca(F));
        % Copy the selected subplot title
        tmp = get(h1,'title'); tmp = tmp.String;
        % Set title to the new figure
        title(gca(F), tmp);
        
        % Set axis to the new figure
        tmp = get(h1,'XLim');
        xlim(gca(F), tmp)
        tmp = get(h1,'YLim');
        ylim(gca(F), tmp)
    case 'alt'
        delete(h1);
end
end


