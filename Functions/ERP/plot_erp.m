function plot_erp(cnfg,ftdata)

%% Plot ERP responses
%
% Syntax:
%    plot_erp(cnfg,ftdata);
%
% Inputs:
%    cfg - Structure of parameters:
%
%       channel  - Channels to plot
%       trigger  - triggers to plot in the same figure
%       latency  - latency to plot
%       dosave   - logical. True/false save/not save the results
%       outpath  - string. Path to save the results if dosave=true
%       plotfig  - logical. Plot the figure with the results.
%       infosave - string to include in the saved filed
%       ampyaxis - 'individual', 'group'. default='group'
%
%   ftdata - Data in format field trip
%
% See also: plot_timelock_sem

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug. 2021; Last revision: 11-Mar-2025

% Change log
% 27/01/26 - The figure when clicking has neg polarity on top of y axis


%% Initialization

if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'plotfig'), cnfg.plotfig=true; end
if ~isfield(cnfg,'latency'), cnfg.latency=ftdata.time{1}([1 end]); end
if ~isfield(cnfg,'channel')
    error('It is mandatory to select channels'), end
if ~isfield(cnfg,'trigger')
    cnfg.trigger = unique(ftdata.trialinfo); end
    %error('It is mandatory to select trigggers'), end
if ~isfield(cnfg,'outpath') && cnfg.dosave
    error('Outpath has not been specified to save the results'), end
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end
if ~isfield(cnfg,'ampyaxis'), cnfg.ampyaxis='group'; end

%% Main workflow

%Creates outpath folder
if cnfg.dosave && ~exist(cnfg.outpath,'dir')
    mkdir(cnfg.outpath)
end

cfg=[];
cfg.channel = cnfg.channel;
ftdata = ft_selectdata(cfg, ftdata);

h=figure;
nrows = ceil(sqrt(length(cfg.channel)));
ncols = ceil(length(cfg.channel)/nrows);

colors = [0         0.4470    0.7410;...
          0.8500    0.3250    0.0980;...
          0.9290    0.6940    0.1250;...
          0.4940    0.1840    0.5560;...
          0.4660    0.6740    0.1880;...
          0.3010    0.7450    0.9330;...
          0.6350    0.0780    0.1840];

max_val = 0;
min_val = 0;
for tr=1:length(cnfg.trigger)
    
    cfg=[];
    cfg.trials  = find(ftdata.trialinfo==cnfg.trigger(tr));
    ftdata_aux = ft_selectdata(cfg, ftdata);
    
    cfg=[];
    cfg.keeptrials='yes';
    timelock = ft_timelockanalysis(cfg, ftdata_aux);
    
    timelock = baselinecorrection(timelock);
    
    resp_max = mean(timelock.trial,1) + ... 
               1.96*std(timelock.trial,[],1)/sqrt(size(timelock.trial,1));
    max_val  = max(max(resp_max(:)),max_val);
    resp_min = mean(timelock.trial,1) - ...
               1.96*std(timelock.trial,[],1)/sqrt(size(timelock.trial,1));    
    min_val  = min(min(resp_min(:)),min_val);
    for iter=1:length(cnfg.channel)
        G(iter) = subplot(nrows,ncols,iter);
        hold on,
        
        cfg=[];
        cfg.channel = iter;
        color_tr = mod(tr-1,7)+1;
        cfg.color   = colors(color_tr,:);
        plot_timelock_sem(cfg,timelock)
        
        %set(h,'UserData',iter);
        %set(h,'HitTest', 'off'); % Disable content selection
        %G(iter).ButtonDownFcn = @newFigure1;
        hold off
    end
end

for iter=1:length(cnfg.channel)
    G(iter) = subplot(nrows,ncols,iter);
    hold on,
    
    title(['ERP: ' ftdata.label{iter}]);
    plot(cnfg.latency,[0 0],'k')
    if strcmp(cnfg.ampyaxis,'group')
        axis([cnfg.latency min_val*1.2 max_val*1.2])
    end
    set(h,'UserData',iter);
    set(h,'HitTest', 'off'); % Disable content selection
    set(gca, 'YDir','reverse')

    G(iter).ButtonDownFcn = @newFigure1;
    hold off
end
% Maximize the figure window to fill the screen
set(h, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

if cnfg.dosave
    savefig(h,[cnfg.outpath 'ERP_resp' cnfg.infosave])
    saveas(h,[cnfg.outpath 'ERP_resp' cnfg.infosave '.png'])
end
if ~cnfg.plotfig
    close(h)
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
        set(gca, 'YDir','reverse')
    case 'alt'
        delete(h1);
end
end

function timelock=baselinecorrection(timelock)
    t1 = find(timelock.time<0,1);
    t2 = find(timelock.time>0,1);
    
    for nch = 1:size(timelock.trial,2)
        timelock.trial(:,nch,:) = timelock.trial(:,nch,:) - mean(timelock.trial(:,nch,t1:t2),3);
    end
end

