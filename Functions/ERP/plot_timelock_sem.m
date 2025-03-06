
function plot_timelock_sem(cnfg,timelock)

% Plot timelock signal of selected channel with confidant interval.
% Also works with data matrix.
%
% USE:
% plot_timelock_sem(cfg,timelock)
%
% INPUT:
% cfg.channel;
% cfg.color: integer (1, 2...) for default colors
%            RGB 
%            Default: gray [0.8 0.8 0.8]       
% cfg.time; % If data in matrix form
% cfg.sig_val; %a binary vector with significant time-points can be used to
%              plot green squares of significance
%

% Log:
% 20.02.2025: Added sig_val
% 28.11.2024: Added default colors
% 03.09.2024: Corrected an error when cfg.time was a column and not a row

% Default Matlab colors
%[0 0.4470 0.7410]	"#0072BD"	
%[0.8500 0.3250 0.0980]	"#D95319"	
%[0.9290 0.6940 0.1250]	"#EDB120"	
%[0.4940 0.1840 0.5560]	"#7E2F8E"	
%[0.4660 0.6740 0.1880]	"#77AC30"	
%[0.3010 0.7450 0.9330]	"#4DBEEE"	
%[0.6350 0.0780 0.1840]

if ~isfield(cnfg,'color')
    color=[0.8 0.8 0.8];
else
    color=cnfg.color;
end
if length(color)==1 % Use default colors
    switch(rem(color,7))
        case 1
            color = [0 0.4470 0.7410];
        case 2
            color = [0.8500 0.3250 0.0980];
        case 3
            color = [0.9290 0.6940 0.1250];
        case 4
            color = [0.4940 0.1840 0.5560];
        case 5
            color = [0.4660 0.6740 0.1880];
        case 6
            color = [0.3010 0.7450 0.9330];
        case 0
            color = [0.6350 0.0780 0.1840];
    end
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
%curve1 = MU + SIGMA;
%curve2 = MU - SIGMA;

ttime = [timey, fliplr(timey)];
inBetween = [curve1, fliplr(curve2)];

hold on
P = fill(ttime, inBetween, color,'LineStyle','none');
alpha(P,0.5)
plot(timey,MU,'Color',color,'Linewidth',1);
hold off        

%%% Plot significance
if isfield(cnfg,'sig_val')
    if length(cnfg.sig_val)~=length(timey)
        error('To plot significance the sig_val vector must have the same length as the ''time'' vector')
    end

    starts=1+find(diff(cnfg.sig_val)==1);
    stops=1+find(diff(cnfg.sig_val)==-1);
    if ~isempty(starts)
        if(stops(end)<starts(end))
            starts=starts(1:end-1);
        end
        if(starts(1)>stops(1))
            stops=stops(2:end);
        end

        %Display the significant values
        hold on
        ax=axis;
        for k=1:length(starts)
            p=patch(timey([starts(k) starts(k) stops(k) stops(k)]), [ax(3) ax(4) ax(4) ax(3)],'g','LineStyle','none');
            set(p,'FaceAlpha',0.2, 'HitTest', 'off')
        end
    end
end





        