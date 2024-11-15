function [pred,stats_TRF]=ft_trf_forward(cnfg, ftdata, feat)

% Function to compute TRF using mTRF
%
% USE:
%   ft_trf_forward(cfg, ftdata, feat)
%
% INPUT:
%   cfg - struct with parameters:
%   Ntr     - If continuous data, create trials (Def=5).
%   tmin    - time prestim in ms.
%   tmax    - time poststim in ms.
%   lambda  - lambda to test (Def=10.^(-15:0.5:15));
%   doplot      - true/false (Def: true)
%   dosave      - true/false (Def: false)
%   outpath     - string with path 
%   infosave    - string with filename

% OUTPUT:
%
% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Mar. 2024; Last revision: 12-Nov-2024

%Change log
% 

% To do:
% check that ftdata and feat have the same size

%% Parameters

if ~isfield(cnfg,'Ntr'), cnfg.Ntr   = 5; end
if ~isfield(cnfg,'tmin'), error('Parameter tmin is mandatory'); end
if ~isfield(cnfg,'tmax'), error('Parameter tmax is mandatory'); end
if ~isfield(cnfg,'lambda'), cnfg.lambda=10.^(-15:0.5:15); end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if isfield(cnfg,'outpath')
    if ~strcmp(cnfg.outpath(end),'\')
        cnfg.outpath = [cnfg.outpath '\'];
    end
end
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end

if cnfg.dosave
    if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
    if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
end

%% Prepare data

fs = ftdata.fsample;
Nch = size(ftdata.trial{1},1);
Nfeat = size(feat,1);

if length(ftdata.trial)==1 % No trials, I must create them
    Ntr = cnfg.Ntr;
    stim = cell(1,Ntr);
    resp = cell(1,Ntr);
    eeg = ftdata.trial{1};
    nsamples = floor(size(eeg,2)/Ntr);
    for tri=1:Ntr
        stim{tri} = feat(:,(tri-1)*nsamples+1:tri*nsamples)';
        resp{tri} = eeg(:,(tri-1)*nsamples+1:tri*nsamples)';
    end
else % data already in trials
    Ntr = length(ftdata.trial);
    stim = cell(1,Ntr);
    resp = cell(1,Ntr);
    for tri=1:Ntr
        stim{tri} = feat{tri}';
        resp{tri} = ftdata.trial{tri}';
    end
end

resp_freq = cell(Nch,Ntr);
for chi=1:Nch
    for tri=1:Ntr
        resp_freq{chi,tri} = resp{tri}(:,chi)';
    end
end

%% TRF - Crossvalidation
% I find the best lambda (ridge parameter) on the whole dataset. 
% Then I train the model on all trials expect one and test in the last trial

% Compute response models
% [stats,t] = mTRFcrossval(stim,resp,fs,Dir,tmin,tmax,lambda,varargin)
% model = mTRFtrain(stim,resp,fs,Dir,tmin,tmax,lambda,varargin)
Dir = 1; %forward model
tmin = cnfg.tmin; %in ms
tmax = cnfg.tmax; %in ms
lambda = cnfg.lambda;

for chi=1:Nch
    [stats{chi},t{chi}] = mTRFcrossval(stim,resp_freq(chi,:),fs,Dir,tmin,tmax,lambda);
    [~,p]=max(mean(mean(stats{chi}.r),3));
    lambda_model(chi) = lambda(p);
end

%figure, hold on,
%plot(log(lambda),zscore(mean(mean(stats{chi}.r),3)))
%plot(log(lambda),zscore(mean(mean(stats{chi}.err),3)))

%% Train + Test

model = cell(Nch,Ntr);
pred = cell(Nch,Ntr);
stats_TRF = cell(Nch,Ntr);
for chi = 1:Nch
    for tri = 1:Ntr
        data_train = 1:Ntr;
        data_train(tri) = []; %This trial is for test
        data_test = tri;
        model{chi,tri} = mTRFtrain(stim(data_train),resp_freq(chi,data_train),fs,Dir,tmin,tmax,lambda_model(chi));
        [pred{chi,tri},stats_TRF{chi,tri}] = mTRFpredict(stim(data_test),resp_freq(chi,data_test),model{chi,tri});
    end
end

%% Plot and Save data

% rsquare = zeros(Nch+1,1);
% for chi = 1:Nch
%     for tri = 1:Ntr
%         rsquare(chi,1:end-1) = rsquare(chi,1:end-1) + stats_TRF{chi,tri}.r.^2;
%     end
% end
% rsquare = rsquare/Ntr;
% 
% freq_axis = cnfg.freq(1:end-1) + diff(cnfg.freq)/2;
% ylabels = cell(1,length(freq_axis));
% for fi=1:length(freq_axis)
%     ylabels{fi} = num2str(round(freq_axis(fi),1));
% end
% for fi=1:length(cnfg.freq)
%     ylabels2{fi} = num2str(round(cnfg.freq(fi),1));
% end
% 
% % h1=figure; hold on,
% % pcolor(1:Nch+1,cnfg.freq,rsquare')
% % yticks(round(cnfg.freq,1))
% % %yticklabels(ylabels2)
% % xticks(1.5:1:Nch+0.5)
% % xticklabels(ftdata.label)
% % axis([1 Nch+1 cnfg.freq(1) cnfg.freq(end)])
% % ylabel('Frequency (Hz)')
% % colorbar
% % colormap("jet")
% % title('TRF Speech-envelope --> Brain')
% 
% h2=figure; hold on,
% pcolor(1:Nch+1,log(cnfg.freq),rsquare')
% yticks(log(cnfg.freq))
% yticklabels(ylabels2)
% xticks(1.5:1:Nch+0.5)
% xticklabels(ftdata.label)
% axis([1 Nch+1 log(cnfg.freq(1)) log(cnfg.freq(end))])
% ylabel('Frequency (Hz)')
% colorbar
% colormap("jet")
% title('TRF Speech-envelope --> Brain')
% 
% if cnfg.dosave
%     %savefig(h1,[cnfg.outpath 'TRF_freq' cnfg.infosave ])
%     %saveas(h1,[cnfg.outpath 'TRF_freq' cnfg.infosave '.png'])
%     savefig(h2,[cnfg.outpath 'TRF_freq_log' cnfg.infosave ])
%     saveas(h2,[cnfg.outpath 'TRF_freq_log' cnfg.infosave '.png'])
% end
% if ~cnfg.doplot, close(h1), close(h2), end
% % Save data
% if cnfg.dosave
%     save([cnfg.outpath 'TRF_freq' cnfg.infosave],'model','stats_TRF','stats')
% end
% 
% %%% Plot TRF model
% 
% % TRF_model = cell(1,Nch);
% % time_model = model{1}.t;
% % Npt = length(time_model);
% % for chi=1:Nch
% %     TRF_model{chi} = zeros(Nfeat,Npt,length(cnfg.freq)-1);
% %     for tri=1:Ntr
% %         TRF_model{chi} = TRF_model{chi} + model{chi}.w;
% %     end
% %     TRF_model{chi} = TRF_model{chi}/Ntr;
% % end
% % 
% % for feati=1:Nfeat
% %     hfeat = figure;
% %     nrows = ceil(sqrt(Nch));
% %     ncols = ceil(Nch/nrows);
% %     for chi=1:Nch
% %         subplot(nrows,ncols,chi), hold on
% %         for fi=1:length(cnfg.freq)-1
% %             plot(time_model,squeeze(TRF_model{chi}(feati,:,fi)),...
% %                 'LineWidth',1,...
% %                 'DisplayName',['Freq: ' num2str(round(cnfg.freq(fi),1)) '-'...
% %                 num2str(round(cnfg.freq(fi+1),1)) ' Hz'])
% %         end
% %     legend
% %     plot(time_model([1 end]),[0 0],'k')
% %     title(ftdata.label{chi})
% %     end
% % end

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
            XLim = get(h1,'XLim');
            YLim = get(h1,'YLim');
            set(gca, 'XLim',XLim)
            set(gca, 'YLim',YLim)
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