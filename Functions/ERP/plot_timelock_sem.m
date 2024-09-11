
function plot_timelock_sem(cnfg,timelock)

% Plot timelock signal of selected channel with confidant interval.
% Also works with data matrix.
%
% USE:
% plot_timelock_sem(cfg,timelock)
%
% INPUT:
% cfg.channel;
% cfg.color;
% cfg.time; % If data in matrix form
%

% Log:
% 03.09.2024: Corrected an error when cfg.time was a column and not a row

if ~isfield(cnfg,'color')
    color=[0.8 0.8 0.8];
else
    color=cnfg.color;
end

if isstruct(timelock)
    data=squeeze(timelock.trial(:,cnfg.channel,:));
    timey = timelock.time; 
else
    data=timelock;
    timey=cnfg.time(:)';
end

MU = mean(data);
SIGMA = std(data);

curve1 = MU + 1.96*SIGMA/sqrt(size(data,1));
curve2 = MU - 1.96*SIGMA/sqrt(size(data,1));
%curve1 = MU + SIGMA/sqrt(size(data,1));
%curve2 = MU - SIGMA/sqrt(size(data,1));

ttime = [timey, fliplr(timey)];
inBetween = [curve1, fliplr(curve2)];

hold on
P = fill(ttime, inBetween, color,'LineStyle','none');
alpha(P,0.5)
plot(timey,MU,'Color',color,'Linewidth',1);
hold off        
        
        