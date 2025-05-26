function snr_erp(cnfg,ftdata)

%% Compute SNR of ERPs
% It is computed as SNR= erp/std(baseline)

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Mar. 2025; Last revision: 17-Mar-2025


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

for tr=1:length(cnfg.trigger)    
    cfg=[];
    cfg.trials  = find(ftdata.trialinfo==cnfg.trigger(tr));
    ftdata_aux = ft_selectdata(cfg, ftdata);
    
    cfg=[];
    cfg.keeptrials = 'yes';
    timelock = ft_timelockanalysis(cfg, ftdata_aux);
    timelock=baselinecorrection(timelock);
    timelock.var = squeeze(std(timelock.trial,[],1));
    timelock.avg = squeeze(mean(timelock.trial,1));
    
    tini=find(timelock.time<=0,1);
    tend=find(timelock.time>=0,1); tend=tend-1;
    
    N=mean(timelock.var(:,tini:tend),2);
    N=repmat(N,1,size(timelock.var,2));
    timelock.avg = abs(timelock.avg) ./ N;
    
    for iter=1:length(cnfg.channel)
        G(iter) = subplot(nrows,ncols,iter);
        hold on,
        plot(timelock.time,timelock.avg(iter,:))
        
        set(h,'UserData',iter);
        set(h,'HitTest', 'off'); % Disable content selection
        G(iter).ButtonDownFcn = @newFigure1;
        title(['ERP: ' ftdata.label{iter}]);
        hold off
    end
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
