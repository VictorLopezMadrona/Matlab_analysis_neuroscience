function evok_resp_cmp=cmp_trig(cnfg,ftdata)

%% Compare the ERP of two triggers
% 
% Syntax:  
%    evok_resp=cmp_trig(cfg,ftdata);
%
% Inputs:
%    cfg - Structure of parameters:
% 
%       eventvalue        = vector with the two triggers to compare.
%
%       channel           = select only some channels. Optional
%
%       minsamples        = integer, the minimum number of consecutive 
%                           samples allowed to be detected as significant. 
%                           Default is 2, where segments of at least two 
%                           consecutive samples will be detected.
%
%       stats             = string. Method to estimate significance.
%                           Options:
%                           'lfdr' (default), 'perm'.
%
%       Nsurro            = number of surrogates for permutation
%                           By default is 100.
%
%       dosave            = logical. True/false save/not save the results
%
%       infosave          = string to include in the saved file
%
%       plotfig           = logical. Plot the results
%
%       outpath           = string. Path to save the results if dosave=true
%
%   ftdata - Data in format field trip
%
% Outputs:
%   evok_resp_cmp
%
% See also:

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Nov. 2019; Last revision: 9-Dec-2020

% Change_log:
% 9-Dec-2020: included 'infosave'

%% Initialization
if ~isfield(cnfg,'minsamples'), cnfg.minsamples=2; end
if ~isfield(cnfg,'latency'), cnfg.latency=[0 ftdata.time{1}(end)]; end
if ~isfield(cnfg,'stats'), cnfg.stats='lfdr'; end 
if ~isfield(cnfg,'Nsurro'), cnfg.Nsurro=100; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end
if ~isfield(cnfg,'plotfig'), cnfg.plotfig=true; end
if ~isfield(cnfg,'outpath') && cnfg.dosave
    error('Outpath has not been specified to save the results')
end

if ~isfield(cnfg,'eventvalue')
    error('Param "eventvalue" must be defined')
end
if length(cnfg.eventvalue)~=2
    error('Comparison must be done between two conditions')
end

%% Main workflow

if isfield(cnfg,'channel')
    cfg=[];
    cfg.channel = cnfg.channel;
    ftdata=ft_selectdata(cfg,ftdata);
else
    cnfg.channel = 'all';
end

eventvalue = cnfg.eventvalue;

% Keep only trials to analyze and remove other triggers
rmv_trl=find(ftdata.trialinfo~=eventvalue(1) & ftdata.trialinfo~=eventvalue(2));
ftdata.time(rmv_trl)=[];
ftdata.trial(rmv_trl)=[];
if isfield(ftdata,'sampleinfo')
    ftdata.sampleinfo(rmv_trl,:)=[]; end
ftdata.trialinfo(rmv_trl)=[];

cfg=[];
cfg.eventvalue = cnfg.eventvalue;
data2Chr = ftdata2evxic(cfg,ftdata);
    
cfg=[];
cfg.latency = [0 ftdata.time{1}(end)];
ftdata_aux = ft_selectdata(cfg,ftdata);
    
cfg=[];
cfg.eventvalue=cnfg.eventvalue;
if strcmp(cnfg.stats,'perm')
    Nsurro=cnfg.Nsurro;
    disp(['Computing ' num2str(Nsurro) ' surrogates...'])
    Nch=length(ftdata.label);
    Nsamples=zeros(Nch,Nsurro);
    
    %%% Pixel-based correction
    SURRO_th = zeros(Nsurro, 2, Nch);
    for s=1:Nsurro      
        s
        ftdata_aux.trialinfo=ftdata_aux.trialinfo(randperm(length(ftdata_aux.trialinfo)));
        data_surro=ftdata2evxic(cfg,ftdata_aux);
        for ch = 1:Nch
            data1 = data_surro{1,ch};
            data2 = data_surro{2,ch};
            
            data_diff = mean(data1) - mean(data2);
            SURRO_th(s,1,ch) = max(data_diff);
            SURRO_th(s,2,ch) = min(data_diff);
        end
        %Nsamples(:,s) = aw_signif_comps_clus(data_surro);
    end
    
    cfg=[];
    cfg.latency = [0 ftdata.time{1}(end)];
    ftdata_aux = ft_selectdata(cfg,ftdata);
    cfg=[];
    cfg.eventvalue=cnfg.eventvalue;
    data_surro=ftdata2evxic(cfg,ftdata_aux);
    sign_pixel = zeros(Nch,size(data_surro{1},2));
    for ch = 1:Nch
        data1 = data_surro{1,ch};
        data2 = data_surro{2,ch};
        data_diff = mean(data1) - mean(data2);
        MU_max = mean(squeeze(SURRO_th(:,1,ch)));
        MU_min = mean(squeeze(SURRO_th(:,2,ch)));
        SIGMA_max = std(squeeze(SURRO_th(:,1,ch)));
        SIGMA_min = std(squeeze(SURRO_th(:,2,ch)));
        
        xmax = norminv(0.975,MU_max,SIGMA_max);
        xmin = norminv(0.025,MU_min,SIGMA_min);
        
        sign_pixel(ch,:) = (data_diff>=xmax) + (data_diff<=xmin);
    end
        
%     %For each surrogate, keep the maximum size independently of the channel
%     minsamples=ceil(norminv(0.95,mean(max(Nsamples)),std(max(Nsamples)))); %Without transforming the data
    
    %For each channel, keep the maximum between surrogates and then the
    %maximum between channels (less restrictive)
    for chi=1:Nch
        minsamples(chi)=ceil(norminv(0.95,mean(Nsamples(chi,:)),std(Nsamples(chi,:)))); %Without transforming the data
    end
    
    minsamples=max(minsamples);
    disp(['Minimum number of consecutive samples = ' num2str(minsamples)])
else
    minsamples=cnfg.minsamples;
end

cfg=[];
cfg.doplot     = 0;
cfg.minsamples = minsamples;
cfg.latency    = cnfg.latency;
cfg.time       = ftdata.time;
cfg.stats      = cnfg.stats;
[signif_comp_idx, tt, t_signif]=aw_signif_comps2(cfg, data2Chr);


%Creates outpath folder
if cnfg.dosave && ~exist(cnfg.outpath,'dir')
     mkdir(cnfg.outpath)
end

if isempty(signif_comp_idx)
    warning('No significant responses were found')
else
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    plotCompStatsNAverage_inter(data2Chr(1,signif_comp_idx), ftdata.label(signif_comp_idx), ftdata.time{1}, [], t_signif(signif_comp_idx,:), [0 114 189]/255);
    plotCompStatsNAverage_inter(data2Chr(2,signif_comp_idx), ftdata.label(signif_comp_idx), ftdata.time{1}, [], t_signif(signif_comp_idx,:), [217 83 25]/255);
    if cnfg.dosave
    savefig(h,[cnfg.outpath 'EvokResp_comp_trig' num2str(eventvalue(1)) '-' num2str(eventvalue(2)) cnfg.infosave])
    saveas(h,[cnfg.outpath 'EvokResp_comp_trig' num2str(eventvalue(1)) '-' num2str(eventvalue(2)) cnfg.infosave '.png'])
    end
    if ~cnfg.plotfig, close(h), end
end
evok_resp_cmp.ICsig = signif_comp_idx';
evok_resp_cmp.t_signif = t_signif;
evok_resp_cmp.tt = tt;
evok_resp_cmp.time = ftdata.time{1};
evok_resp_cmp.label = ftdata.label;
evok_resp_cmp.channel = cnfg.channel;
evok_resp_cmp.cnfg = cnfg;
if cnfg.dosave
    save([cnfg.outpath 'Evok_resp_cmp_trig' num2str(eventvalue(1)) '-' num2str(eventvalue(2)) cnfg.infosave],'evok_resp_cmp')
end



