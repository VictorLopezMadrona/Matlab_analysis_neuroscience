
function [hsample,psample,pval_corrected]=compare2signal(cnfg,data,data2)

%% Plot ERP responses
%
% Syntax:
%    compare2signal(cnfg,data,data2);
%
% Inputs:
%    cfg - Structure of parameters:
%
%       color    - RGB values with the color of the plot (not working)
%       dosave   - logical. True/false save/not save the results (not working)
%       outpath  - string. Path to save the results if dosave=true (not working)
%       plotfig  - logical. Plot the resultant figure. (not working)
%       infosave - string to include in the saved filed (not working)
%       mode     - type of stats (1, 2, 3...) (default FDR on each channel: 3)
%       pval     - Def 0.05
%       Nperm    - Def 100
%
%   data - matrix with the data [Nsamples,Nch,Nrep] or [Nsamples,Nrep]
%
% MODE 1: FDR on all the pvalues
% MODE 2: FDR on each channel - we keep pID
% MODE 3: FDR on each channel - we keep pN (more restrictive)
% MODE 4: Permutations uncorrected
% MODE 5: Permutations pixel-based correction
%
% Outputs:
%       hsample - 0 non-significative; 1 significative (corrected)
%       psample - pvalue of each sample
%       pval_corrected - corrected p-value
%
% See also:

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Feb. 2024; Last revision: 2-Feb-2024

% Change-log:
% 21/08/2024: Included correction using permutations

%% Initialization

if ~isfield(cnfg,'infosave'), cnfg.infosave = ''; end
if ~isfield(cnfg,'dosave'),   cnfg.dosave   = false; end
if ~isfield(cnfg,'plotfig'),  cnfg.plotfig  = true; end
if ~isfield(cnfg,'color'),    cnfg.color    = [0 114 189]/255; end
if ~isfield(cnfg,'outpath') && cnfg.dosave
    error('Outpath has not been specified to save the results'), end
if ~isfield(cnfg,'pval'), cnfg.pval = 0.05; end

%%% MODE 1: FDR on all the pvalues
%%% MODE 2: FDR on each channel - we keep pID
%%% MODE 3: FDR on each channel - we keep pN (more restrictive)
if ~isfield(cnfg,'mode'), cnfg.mode = 3; end
if ~isfield(cnfg,'Nperm'), cnfg.Nperm = 100; end


%% Initial testing

%Creates outpath folder
if cnfg.dosave && ~exist(cnfg.outpath,'dir')
    mkdir(cnfg.outpath)
end

if length(size(data))==2 %Only one channel
    [Nsamples1,~] = size(data);
    Nch1 = 1;
else
    [Nsamples1,Nch1,~] = size(data);
end

if length(size(data2))==2 %Only one channel
    [Nsamples2,~] = size(data2);
    Nch1 = 2;
else
    [Nsamples2,Nch2,~] = size(data2);
end

if Nsamples1~=Nsamples2
    error('The number of samples in both datasets must be the same')
end
if Nch1~=Nch2
    error('The number of channels in both datasets must be the same')
end

%% FDR

if cnfg.mode==1 || cnfg.mode==2 || cnfg.mode==3
pval=zeros(Nsamples1,Nch1);
tval=zeros(Nsamples1,Nch1);

for chi=1:Nch1
    [~,pval(:,chi),~,stats]=ttest2(squeeze(data(:,chi,:))',squeeze(data2(:,chi,:))');
    tval(:,chi)=stats.tstat;
end

if cnfg.mode==1
    [pID,pN] = fdr(pval,cnfg.pval);
end
if cnfg.mode==2 || cnfg.mode==3
    for chi=1:Nch1
        [pID_aux,pN_aux]= fdr(pval(:,chi),cnfg.pval);
        if isempty(pID_aux), pID(chi)=0; else, pID(chi)=pID_aux; end
        if isempty(pN_aux), pN(chi)=0; else, pN(chi)=pN_aux; end
    end
end

psample = pval;
if cnfg.mode==1
    pval_corrected = pID;
    hsample = pval<=pval_corrected;
elseif cnfg.mode==2
    pval_corrected = pID;
    for chi=1:Nch1
        hsample(:,chi) = pval(:,chi)<=pval_corrected(chi);
    end
elseif cnfg.mode==3
    pval_corrected = pN;
    for chi=1:Nch1
        hsample(:,chi) = pval(:,chi)<=pval_corrected(chi);
    end
end
end

%% PERMUTATIONS
if cnfg.mode==4 || cnfg.mode==5

    pval=zeros(Nsamples1,Nch1);
    tval=zeros(Nsamples1,Nch1);
    tval_perm=zeros(Nsamples1,Nch1,cnfg.Nperm);
    [Nsamples,Nch1,Ntr1] = size(data);
    [~,~,Ntr2] = size(data);
    diff_perm = zeros(Nsamples,Nch1,cnfg.Nperm);
    data_all  = zeros(Nsamples,Nch1,Ntr1+Ntr2);
    data_all(:,:,1:Ntr1) = data;
    data_all(:,:,Ntr1+1:end) = data2;

    for chi=1:Nch1
        [~,pval(:,chi),~,stats]=ttest2(squeeze(data(:,chi,:))',squeeze(data2(:,chi,:))');
        tval(:,chi)=stats.tstat;
    end
    
    %data_diff = mean(data,3)-mean(data2,3);
    for permi = 1:cnfg.Nperm
        rp = randperm(Ntr1+Ntr2);
        data1_aux = data_all(:,:,rp(1:Ntr1));
        data2_aux = data_all(:,:,rp(Ntr1+1:end));
        for chi=1:Nch1
            [~,~,~,stats]=ttest2(squeeze(data1_aux(:,chi,:))',squeeze(data2_aux(:,chi,:))');
            tval_perm(:,chi,permi)=stats.tstat;
        end
        %diff_perm(:,:,permi) = mean(data1_aux,3)-mean(data2_aux,3);
    end
    
    if cnfg.mode==4 
    %%% Uncorrected
    MU_perm = mean(tval_perm,3);
    SIGMA_perm = std(tval_perm,0,3);
    psample = normcdf(tval,MU_perm,SIGMA_perm);
    hsample = (psample>=0.975) + (psample<=0.025);
    
    elseif cnfg.mode==5
    %%% Pixel-based correction
    diff_perm_min = min(tval_perm);
    diff_perm_max = max(tval_perm);

    for chi=1:Nch1
        MU_perm = mean(squeeze(diff_perm_min(:,chi,:)));
        SIGMA_perm = std(squeeze(diff_perm_min(:,chi,:)));
        psample_min(:,chi) = normcdf(squeeze(tval(:,chi)),MU_perm,SIGMA_perm);
        
        MU_perm = mean(squeeze(diff_perm_max(:,chi,:)));
        SIGMA_perm = std(squeeze(diff_perm_max(:,chi,:)));
        psample_max(:,chi) = normcdf(squeeze(tval(:,chi)),MU_perm,SIGMA_perm);
    
        hsample = (psample_max>=0.975) + (psample_min<=0.025);
        psample = min(1-psample_max,psample_min);
    end
    end

    pval_corrected = NaN;
end


% if nargin==3
%     %Select some channels
%     aw_data2.trial = aw_data2.trial(cnfg.channel);
%     aw_data2.label = aw_data2.label(cnfg.channel);
%     %Change format of aw_data.trial
%     data = aw_data2.trial;
%     aw_data2.trial = [];
%     Nch = length(data);
%     Ntr = size(data{1},1);
%     for chi=1:Nch
%         for tri=1:Ntr
%             aw_data2.trial(tri,chi,:) = data{chi}(tri,:);
%         end
%     end
% end
% 
% if isfield(cnfg,'pval')
%     if nargin==2, sig_time = compute_significance(cnfg.pval,aw_data); end
%     if nargin==3, sig_time = compute_significance(cnfg.pval,aw_data,aw_data2); end
% end
% 
% h=figure;
% nrows = ceil(sqrt(length(cnfg.channel)));
% ncols = ceil(length(cnfg.channel)/nrows);
% 
% max_val = 0;
% min_val = 0;
% 
% %aw_data = baselinecorrection(aw_data);
% 
% resp_max = mean(aw_data.trial,1) + ...
%     1.96*std(aw_data.trial,[],1)/sqrt(size(aw_data.trial,1));
% max_val  = max(max(resp_max(:)),max_val);
% resp_min = mean(aw_data.trial,1) - ...
%     1.96*std(aw_data.trial,[],1)/sqrt(size(aw_data.trial,1));
% min_val  = min(min(resp_min(:)),min_val);
% for iter=1:length(cnfg.channel)
%     G(iter) = subplot(nrows,ncols,iter);
%     hold on,
%     
%     cfg=[];
%     cfg.channel = iter;
%     cfg.color   = cnfg.color;
%     plot_timelock_sem(cfg,aw_data)
%     
%     %set(h,'UserData',iter);
%     %set(h,'HitTest', 'off'); % Disable content selection
%     %G(iter).ButtonDownFcn = @newFigure1;
%     hold off
% end
% 
% if nargin==3
%     %aw_data2 = baselinecorrection(aw_data2);
%     resp_max = mean(aw_data2.trial,1) + ...
%         1.96*std(aw_data2.trial,[],1)/sqrt(size(aw_data2.trial,1));
%     max_val  = max(max(resp_max(:)),max_val);
%     resp_min = mean(aw_data2.trial,1) - ...
%         1.96*std(aw_data2.trial,[],1)/sqrt(size(aw_data2.trial,1));
%     min_val  = min(min(resp_min(:)),min_val);
%     for iter=1:length(cnfg.channel)
%         G(iter) = subplot(nrows,ncols,iter);
%         hold on,
%         cfg=[];
%         cfg.channel = iter;
%         cfg.color   = [217 83 25]/255;
%         plot_timelock_sem(cfg,aw_data2)
%         hold off
%     end
% end
% 
% for iter=1:length(cnfg.channel)
%     G(iter) = subplot(nrows,ncols,iter);
%     hold on,
%     
%     title(['ERP ' aw_data.label{iter}]);
%     plot(cnfg.latency,[0 0],'k')
%     %axis([cnfg.latency min_val*1.2 max_val*1.2])
%     xlim([cnfg.latency])
%     set(h,'UserData',iter);
%     set(h,'HitTest', 'off'); % Disable content selection
%     set(gca, 'YDir','reverse')
%     
%     G(iter).ButtonDownFcn = @newFigure1;
%     hold off
% end
% 
% % Plot significance
% % Computations to display the statistically significant segments
% for iter=1:length(cnfg.channel)
% starts=1+find(diff([0 sig_time(iter,:) 0])==1);
% stops=1+find(diff([0 sig_time(iter,:) 0])==-1);
% starts=starts-2; stops=stops-2;
% starts(starts<=0)=1;
% stops(stops<=0)=1;
% if ~isempty(starts)
%         if(stops(end)<starts(end))
%             starts=starts(1:end-1);
%         end
%         if(starts(1)>stops(1))
%             stops=stops(2:end);
%         end
%     %Display the significant values for each channels
%     G(iter) = subplot(nrows,ncols,iter);
%     hold on,
%     ax=axis;
%     for k=1:length(starts)
%         p=patch(aw_data.time([starts(k) starts(k) stops(k) stops(k)]), [ax(3) ax(4) ax(4) ax(3)],'g','LineStyle','none');
%         set(p,'FaceAlpha',0.2, 'HitTest', 'off')
%     end
% end
% end
%     
% 
% if cnfg.dosave
%     savefig(h,[cnfg.outpath 'ERP_resp' cnfg.infosave])
%     saveas(h,[cnfg.outpath 'ERP_resp' cnfg.infosave '.png'])
% end
% if ~cnfg.plotfig
%     close(h)
% end
% 
% 
% end
% 
% function newFigure1(h1,~)
% %% Function to act on subplot click
% %         Mouse click: Plots the selected subplot to a new figure
% %  Ctrl + Mouse click: Delete subplot
%     switch get(gcf,'SelectionType')
%         case 'normal'
%             F = figure();
%             copyobj(h1.Children,gca(F));
%             % Copy the selected subplot title
%             tmp = get(h1,'title'); tmp = tmp.String;
%             % Set title to the new figure
%             title(gca(F), tmp);
%             set(gca, 'YDir','reverse')
%             XLim = get(h1,'XLim');
%             YLim = get(h1,'YLim');
%             set(gca, 'XLim',XLim)
%             set(gca, 'YLim',YLim)
%         case 'alt'
%             delete(h1);
%             
%             if false
%                 % Get the list of remaining channels (stored as titles)
%                 fig=flipud(findobj(gcf,'type','axes'));
%                 % Store titles to a cell array
%                 subplot_titles = cell(length(fig),1);
%                 for i=1:length(fig)
%                     tmp = get(fig(i),'title');
%                     subplot_titles(i) = tmp.String;
%                 end
%                 % Coma delimited channels
%                 %sprintf('%s,' , subplot_titles{:})
%             end
%     end
% end
% 
% function sig_samples = compute_significance(pval,data,data2)
% 
% nch = size(aw_data.trial,2);
% p = ones(nch,size(aw_data.time,2));
% for chi=1:nch
%     if nargin==2
%         data_aux = squeeze(aw_data.trial(:,chi,:));
%         [~,p(chi,:)]=ttest(data_aux);
%     else
%         data_aux = squeeze(aw_data.trial(:,chi,:));
%         data_aux2 = squeeze(aw_data2.trial(:,chi,:));
%         [~,p(chi,:)]=ttest2(data_aux,data_aux2);
%     end
% end
% t0=find(aw_data.time>=0,1);
% p(:,1:t0)=1;
% 
% sig_time = zeros(nch,size(aw_data.time,2));
% sig_time(p<=pval)=1;
% end
% 
% % function timelock=baselinecorrection(timelock)
% %     t1 = find(timelock.time<0,1);
% %     t2 = find(timelock.time>0,1);
% %     
% %     for nch = 1:size(timelock.trial,2)
% %         baseline = mean(timelock.trial(:,nch,t1:t2));
% %         baseline = mean(baseline(:));
% %         timelock.trial(:,nch,:) = timelock.trial(:,nch,:) - baseline;
% %     end
% % end
% 
