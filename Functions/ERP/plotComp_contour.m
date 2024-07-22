function plotComp_contour(AnyWaveData1, ChannelsName1, time1, response)
%% Subplots with all trials in an image
% It reorders the trials based on the delay
%
% response(:,1) = order;
% response(:,2) = delay (s);

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Nov. 2019; Last revision: 10-Feb-2021

% Change_log:
% 09-Aug-2021: reorder trials based on response and plot it
% 10-Feb-2021: reorder trials based on delay

% Fix: When no data are given, exit without plotting.
if isempty(AnyWaveData1), return; end

if nargin==4
    flag_resp = 1;
else
    flag_resp = 0;
end

% Extract data shape info
if size(AnyWaveData1,2)>1, AnyWaveData1=AnyWaveData1'; end
nb_comps    = size(AnyWaveData1,1);
%nb_comps    = length(AnyWaveData1);
nb_trials  = size(AnyWaveData1{1},1);

% Compute the number of rows and cols for the subplots
nrows       = ceil(sqrt(nb_comps));
ncols       = ceil(nb_comps/nrows);

for i = 1:nb_comps
    G(i) = subplot(nrows,ncols,i);
       
    %implot=imgaussfilt(AnyWaveData1{i},1);
    implot=AnyWaveData1{i};
    if flag_resp
        implot=implot(response(:,3),:);
    end
    % Reorder delay
    % implot=reorder_trials(implot);
    
    if true % Normalization and smoothing
        implot=zscore(implot,0,2);
        h=[0.05 0.15 0.6 0.15 0.05]';
        implot=imfilter(implot,h);
    end
        
    P=imagesc(time1,1:size(AnyWaveData1{i},1),implot,[-max(abs(implot(:))) max(abs(implot(:)))]);
    colormap(jet)
    
    set(P,'UserData',i);
    set(P, 'HitTest', 'off'); % Disable content selection
    
    G(i).ButtonDownFcn = @newFigure1; % Set click action
    hold on
    plot([0 0],[0 size(AnyWaveData1{i},1)],'k','LineWidth',1.5);
    if flag_resp
        plot(response(response(:,3),2)/1000,1:size(response,1),'k','LineWidth',1.5)
    end
    hold off;
    title(ChannelsName1(i,1));
    xlabel('Time (s)')
    ylabel('# Trial')    
    
end

clear G;

end

function implot=reorder_trials(implot)

% Maximum value is positive or negative?
mimplot=mean(implot,1);
[~,p]=max(abs(mimplot));
if mimplot(p)<0 %maximum is negative
    flag=1;
else
    flag=0;
end
if flag %maximum is negative
    implot=-implot; %I inverse here, and later I inverse again
end

[~,p]=max(implot,[],2); %Find sample with maximum value for each trial
%reorder trials
[~,ord]=sort(p);
implot=implot(ord,:);

if flag %maximum is negative
    implot=-implot; %I inverse here, and later I inverse again
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