function evok_resp=find_sign_trig(cnfg,ftdata)

%% Find significant channels responding to a trigger.
% Analyzes all triggers included in ftdata, considering the time window
% from zero to the end of the trial.
% 
% Syntax:  
%    evok_resp=find_sign_trig(cfg,ftdata);
%
% Inputs:
%    cfg - Structure of parameters:
%       
%       minsamples        = integer, the minimum number of consecutive 
%                           samples allowed to be detected as significant. 
%                           Default is 2, where segments of at least two 
%                           consecutive samples will be detected.
% 
%       totalminsamples   = integer, the total minimum number of samples.
%                           Components with less significant samples will be
%                           discarded. Defaults is minsamples.
%
%       stats             = string. Method to estimate significance. The only
%                           method implemented is 'lfdr' (default).
%
%       dosave            = logical. True/false save/not save the results
%
%       outpath           = string. Path to save the results if dosave=true
%
%       plotfig           = logical. Plot the resultant figure.
%
%       delay_resp        = [Ntrials x 2] Order of the trials and delay of
%                           the response to include in the plot_contour
%
%       infosave          = string to include in the saved filed
%
%       plot_raster       = logical. Compute raster plot
% 
%   ftdata - Data in format field trip
%
% Outputs:
%   evok_resp
%
% See also: plotCompStatsNAverage_inter aw_signif_comps2 cmp_trig

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Nov. 2019; Last revision: 13-Dec-2019

% Change_log:
% 19/08/2024: Added p-value as output
% 10/08/2021: Add option to plot the response in the contourplot
%

% To do:
%   - Surrogates to find minimum number of consecutive samples

%% Initialization
if ~isfield(cnfg,'minsamples'), cnfg.minsamples=2; end
if ~isfield(cnfg,'totalminsamples'), cnfg.totalminsamples=cnfg.minsamples; end
if ~isfield(cnfg,'stats'), cnfg.stats='lfdr'; end %so far, only works with lfdr
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'latency'), cnfg.latency=[0 ftdata.time{1}(end)]; end
if ~isfield(cnfg,'plotfig'), cnfg.plotfig=true; end
if ~isfield(cnfg,'outpath') && cnfg.dosave
    error('Outpath has not been specified to save the results')
end
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end

if ~isfield(cnfg,'plot_raster') && isfield(cnfg,'delay_resp') 
    cnfg.plot_raster=true;
elseif ~isfield(cnfg,'plot_raster') && ~isfield(cnfg,'delay_resp')
    cnfg.plot_raster=false;
end
if isfield(cnfg,'delay_resp')
    if size(cnfg.delay_resp,1)~=length(ftdata.trial)
        warning('The number of responses must the same as number of trials')
        cnfg = rmfield(cnfg,'delay_resp');
    end
end    

%% Main workflow

eventvalue = unique(ftdata.trialinfo);

%Creates outpath folder
if cnfg.dosave && ~exist(cnfg.outpath,'dir')
     mkdir(cnfg.outpath)
end

for ev=1:length(eventvalue)
    
    %Create the delay_resp structure for each 
    if isfield(cnfg,'delay_resp')
        %Select only analyzed events
        clear delay_resp
        delay_resp(:,2)=cnfg.delay_resp(ftdata.trialinfo==eventvalue(ev));
        [~,delay_resp(:,3)]=sort(delay_resp(:,2));
        delay_resp(delay_resp(:,3),1) = 1:size(delay_resp,1);
    end
    
    %Find the min_samples by permuting the baseline
    %if strcmp(cnfg.stats_minsamples,'perm')
%     if false
%     Nsurro = 10;
%     for s=1:Nsurro    
%         cfg=[];
%         cfg.eventvalue = eventvalue(ev);
%         data_aw=ftdata2evxic(cfg,ftdata);
%         for nc=1:length(data_aw)
%             rp=randperm(size(data_aw{nc},2));    
%             for nt=1:size(data_aw{nc},1)    
%                 if rp(nc)>1
%                 data_surro{nc}(nt,:) = [data_aw{nc}(nt,rp(nt):end) data_aw{nc}(nt,1:rp(nt)-1)];
%                 end
%             end
%         end
%         Nsamples(:,s) = aw_signif_comps_clus(data_surro);
%     end
%     end
    
    cfg=[];
    cfg.eventvalue = eventvalue(ev);
    data_aw=ftdata2evxic(cfg,ftdata);
    
    cfg=[];
    cfg.minsamples      = cnfg.minsamples;
    cfg.latency         = cnfg.latency;
    cfg.time            = ftdata.time;
    cfg.stats           = cnfg.stats;
    cfg.totalminsamples = cnfg.totalminsamples;
    [signif_comp_idx, tt, t_signif, pval]=aw_signif_comps2(cfg, data_aw);
   
    if isempty(signif_comp_idx)
        warning('No significant responses were found')
    else
        h=figure('units','normalized','outerposition',[0 0 1 1]);
        plotCompStatsNAverage_inter(data_aw(signif_comp_idx), ftdata.label(signif_comp_idx), ftdata.time{1}, [], t_signif(signif_comp_idx,:));
        if cnfg.dosave
            savefig(h,[cnfg.outpath 'EvokResp_trig' num2str(eventvalue(ev)) cnfg.infosave])
            saveas(h,[cnfg.outpath 'EvokResp_trig' num2str(eventvalue(ev))  cnfg.infosave '.png'])
        end
        if ~cnfg.plotfig
            close(h)
        end
        
        %%% Raster Plot
        if cnfg.plot_raster
        h=figure('units','normalized','outerposition',[0 0 1 1]);
        if isfield(cnfg,'delay_resp')
            plotComp_contour(data_aw(signif_comp_idx), ftdata.label(signif_comp_idx), ftdata.time{1},delay_resp);
        else
            plotComp_contour(data_aw(signif_comp_idx), ftdata.label(signif_comp_idx), ftdata.time{1});
        end
        if cnfg.dosave
            savefig(h,[cnfg.outpath 'EvokResp_contour_trig' num2str(eventvalue(ev)) cnfg.infosave])
            saveas(h,[cnfg.outpath 'EvokResp_contour_trig' num2str(eventvalue(ev))  cnfg.infosave '.png'])
        end
        if ~cnfg.plotfig
            close(h)
        end
        end
    end
    
    for c=1:length(data_aw)
        data_mean(c,:) = mean(data_aw{c});
    end
    evok_resp.ICsig{ev} = signif_comp_idx';
    evok_resp.t_signif{ev} = t_signif;
    evok_resp.data_mean{ev} = data_mean;
    evok_resp.tt{ev} = tt;
    evok_resp.pval{ev} = pval;
end
evok_resp.cnfg = cnfg;
if cnfg.dosave
    save([cnfg.outpath 'Evok_resp' num2str(eventvalue(ev))],'evok_resp')
end

% %% Read SEEG data
% if cnfg.sign_SEEG
% fname = cnfg.datasetSEEG;
% eventvalue = cnfg.eventvalueSEEG;
% 
% cfg=[];
% cfg.dataset=fname;
% cfg.ICname=ICname;
% cfg.trialdef.eventtype  = 'Stimulus';
% cfg.trialdef.eventvalue = eventvalue;
% cfg.trialdef.prestim    = prestim; % in seconds
% cfg.trialdef.poststim   = poststim; % in seconds
% cfg = ft_definetrial(cfg);
% 
% cfg.channel    = {'all'};
% ftdata = ft_preprocessing(cfg);
% 
% % Load ICA matrices
% [nch_meg,samples]=size(ftdata.trial{1});
% if isfield(cnfg,'ICnameSEEG')
%     
%     clear unmixing labels
%     load(cnfg.ICnameSEEG)
%     %Anywave option
%     if exist('unmixing','var') && exist('labels','var')
%     %FT option: src_ica
%     elseif exist('src_ica','var')
%         unmixing = src_ica.unmixing;
%         labels = src_ica.topolabel;
%     else
%         error('Do not recognize ICA input information')
%     end
%     
%     [~,ord]=ismember(ftdata.label,labels);
%     unmixing=unmixing(:,ord);
%     [nic,nch_ic]=size(unmixing);
%     if nch_ic~=nch_meg
%         error('Number of MEG channels and IC matrix must agree')
%     end
% else
%     unmixing = eye(nch_meg);
%     nic = nch_meg;
% end
% 
% %Create matrices with the ICA and trials separately
% %eventvalue=[2 4];
% data_trlica=zeros(nic,samples,length(ftdata.trial));
% for trl_i = 1:length(ftdata.trial)
%     data_trlica(:,:,trl_i) = unmixing*ftdata.trial{trl_i};
% end
% 
% %Create labels
% ICAlabel=cell(nic,1);
% for nic_i=1:nic
%     ICAlabel{nic_i,1} = ['ICA SEEG ' num2str(nic_i)];
% end
% 
% %Create the struct and analyze the data (with Christos code)
% 
% %Creates outpath folder
% if (cnfg.savefigs || cnfg.savedata) && ~exist(cnfg.outpath,'dir')
%      mkdir(cnfg.outpath)
% end
% data_aw=cell(nic,1);
% for nic_i=1:nic
%     data_aw{nic_i}=squeeze(data_trlica(nic_i,:,:))';
% end
% cfg=[];
% cfg.minsamples=cnfg.minsamples;
% [signif_comp_idx, tt, t_signif]=aw_signif_comps(cfg, data_aw);
% %signif_comp_idx(1)=11;
% h=figure('units','normalized','outerposition',[0 0 1 1]);
% plotCompStatsNAverage_inter(data_aw(signif_comp_idx), ICAlabel(signif_comp_idx), ftdata.time{1}, [], t_signif(signif_comp_idx,:));
% savefig(h,[cnfg.outpath 'EvokResp_SEEG_trig' num2str(cnfg.eventvalueMEG)])
% saveas(h,[cnfg.outpath 'EvokResp_SEEG_trig' num2str(cnfg.eventvalueMEG) '.png'])
% 
% evok_resp.SEEG.ICsig = signif_comp_idx;
% evok_resp.SEEG.t_signif = t_signif;
% evok_resp.SEEG.tt = tt;
% evok_resp.SEEG.cnfg = cnfg;
% end
% 
% %save([cnfg.outpath 'Evok_resp' num2str(cnfg.eventvalueMEG)],'evok_resp')

