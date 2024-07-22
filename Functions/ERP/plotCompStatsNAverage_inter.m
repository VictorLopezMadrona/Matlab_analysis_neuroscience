function plotCompStatsNAverage_inter(AnyWaveData1, ChannelsName1, time1, MoyenneSurEpoques1, t_signif, color)
%% Subplots with statistics for each component
% Computes statistic based on each component (like channels or ICA
% components) and subplot them. Each subplot can be plotted separately by
% clicking on the desired plot. Statistics (per component) include:
%   - Average
%   - Type1 Error
%   - Statistically significant (average) samples
%
% Usage:
%   plotCompStatsNAverage_inter(AnyWaveData, ChannelsName, time)
%
%   plotCompStatsNAverage_inter(..., Mean, t_signif)
%
% Tips
%   Part or selected components can be used as:
%       plotCompStatsNAverage_inter(Data1(I), Channels(I), time, MeanData(I,:), t_sig(I,:));
%
%

if nargin==5, color=[0.8 0.8 0.8]; end

% Fix: When no data are given, exit without plotting.
if isempty(AnyWaveData1), return; end

% Extract data shape info
if size(AnyWaveData1,2)>1, AnyWaveData1=AnyWaveData1'; end
nb_comps    = size(AnyWaveData1,1);
%nb_comps    = length(AnyWaveData1);
nb_trials  = size(AnyWaveData1{1},1);

% Compute the number of rows and cols for the subplots
nrows       = ceil(sqrt(nb_comps));
ncols       = ceil(nb_comps/nrows);

% flags for (re)computing statistics
f_recomp_mean     =  ~exist('MoyenneSurEpoques1','var') || isempty(MoyenneSurEpoques1);
% f_recomp_t1error  =  ~exist('pointe_erreur_type1','var');
f_recomp_ttest    =  ~exist('t_signif','var') || isempty(t_signif);

for i = 1:nb_comps
    G(i) = subplot(nrows,ncols,i);
        
    % Compute the mean of each component
    if f_recomp_mean
        MoyenneSurEpoques1(i,:) = mean(AnyWaveData1{i,1},1);
        %MoyenneSurEpoques1(i,:) = mean(AnyWaveData1{i},1);
    end
    % Compute the type1 error
    %if f_recomp_t1error
    pointe_erreur_type1 = std(AnyWaveData1{i,1})/sqrt(nb_trials);
    %end
    
    if false
        P = plot(time1,MoyenneSurEpoques1(i,:));
        
        % Define action on subplot click and disble content selection
        set(P,'UserData',i);
        set(P, 'HitTest', 'off'); % Disable content selection
        G(i).ButtonDownFcn = @newFigure1; % Set click action
        
        hold on
        
        plot(time1,MoyenneSurEpoques1(i,:) + 1.96*pointe_erreur_type1,'r');
        plot(time1,MoyenneSurEpoques1(i,:) - 1.96*pointe_erreur_type1,'r');
        hold off;
    else
        curve1 = MoyenneSurEpoques1(i,:) + 1.96*pointe_erreur_type1;
        curve2 = MoyenneSurEpoques1(i,:) - 1.96*pointe_erreur_type1;

        ttime = [time1, fliplr(time1)]; % Use ; instead of ,
        inBetween = [curve1, fliplr(curve2)]; % Use ; instead of ,
        hold on,
        P = fill(ttime, inBetween, color,'LineStyle','none');
        alpha(P,0.5)
        % Define action on subplot click and disble content selection
        set(P,'UserData',i);
        set(P, 'HitTest', 'off'); % Disable content selection
        G(i).ButtonDownFcn = @newFigure1; % Set click action
        %         plot(time1,MoyenneSurEpoques1(i,:),'r');
        hold on
        plot(time1,MoyenneSurEpoques1(i,:),'Color',color,'Linewidth',1);
        hold off;
    end
    title(ChannelsName1(i,1));
    
    % Locate significant segments (t-test) for each component
    if f_recomp_ttest
        %A = AnyWaveData1{i,1};
        [hhh,p,ci,stats] = ttest( AnyWaveData1{i,1} );
        t_val(i,:)= stats.tstat;
        t_signif(i,:) = hhh;
        %ttestsum(i,:) = sum( t_signif(i,:));
    end
    
    % Computations to display the statistically significant segments
    %FR: Calcul pour affichage des rectangle
    starts=1+find(diff(t_signif(i,:))==1);
    stops=1+find(diff(t_signif(i,:))==-1);
    if ~isempty(starts)
        try
            if(stops(end)<starts(end))
                starts=starts(1:end-1);
            end
            if(starts(1)>stops(1))
                stops=stops(2:end);
            end
        catch
            continue
        end
        %Display the significant values for each channels
        hold on
        ax=axis;
        for k=1:length(starts)
            p=patch(time1([starts(k) starts(k) stops(k) stops(k)]), [ax(3) ax(4) ax(4) ax(3)],'g','LineStyle','none');
            set(p,'FaceAlpha',0.2, 'HitTest', 'off')
        end
    end
end

% Update axis to a tight fit
for i=1:nb_comps
    axis(G(i),[min(time1) max(time1) min(MoyenneSurEpoques1(:))-abs(0.2*min(MoyenneSurEpoques1(:))) max(MoyenneSurEpoques1(:))+abs(0.2*max(MoyenneSurEpoques1(:)))]);
end

clear G;

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
            
        case 'alt'
            delete(h1);
            
            if false
                % Get the list of remaining channels (stored as titles)
                fig=flipud(findobj(gcf,'type','axes'));
                % Store titles to a cell array
                subplot_titles = cell(length(fig),1);
                for i=1:length(fig)
                    tmp = get(fig(i),'title');
                    subplot_titles(i) = tmp.String;
                end
                % Coma delimited channels
                %sprintf('%s,' , subplot_titles{:})
            end
    end
end