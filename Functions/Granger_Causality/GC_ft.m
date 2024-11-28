function [Favg,favg]=GC_ft(cnfg,ftdata)

%% Compute Granger Causality between selected channels in ft structure
%
% USE:
%   GC=GC_ft(cfg,ftdata)
%
% PARAMETERS (cfg.): 
%   
%   latency     - [t_start t_end]
%   channel     - If vector: GC between all selected channels. Ex: [1 5 6];
%               - If cell: GC for each cell. Ex: ch{1}=[1 5 6]; ch{2}=[1 6 7 8];
%   window      - window length for GC in time (s) [Def=5]
%   step        - step between windows in time (s) [Def=window]
%   model_order - Model order (p) for GC [Def=15]
%   resample    - Fs to resample data after filtering
%   Nsurro      - Number of surrogates to statistical significance (Def=0)
%   doplot      - true/false (Def: true)
%   dosave      - true/false (Def: false)
%   outpath     - string with path 
%   infosave    - string with filename
%
%% OUTPUT: Modify this
%   comodulogram = vector of cells:
%       Fs: Sampling Frequency (Hz).
%       bins: Number of bins.
%       f_theta: Information about the frequency theta vector.
%       f_gamma: Information about the frequency gamma vector.
%       MI: (N_gamma_pixels, N_theta_pixels) Matrix with the Modulation.
%           Index of each gamma-theta frequency.
%       CFC: (N_gamma_pixels, N_theta_pixels, 2*bins) Matrix with the
%            average of gamma amplitue during a cycle for each 
%            gamma-theta frequency.
%
% See also: computeGC_MVGC

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Jul 23; Last revision: 06-Jul-2023

%% TO DO:
% - WORK ON SURROGATES. Include them as output of computeGC_MVGC

% 04/10/2023: Corrected error in resampling when Fs was not integer

%% Initial parameters:

if ~isfield(cnfg,'latency')
    cnfg.latency=[ftdata.time{1}(1) ftdata.time{1}(end)]; end
if ~isfield(cnfg,'channel'), cnfg.channel='all'; end
if ~isfield(cnfg,'window'), cnfg.window = 5; end
if ~isfield(cnfg,'step'), cnfg.step = cnfg.window; end
if ~isfield(cnfg,'model_order'), cnfg.model_order = 15; end
if ~isfield(cnfg,'stats'), cnfg.stats = 'cluster'; end
if ~isfield(cnfg,'Nsurro'), cnfg.Nsurro = 0; end
if ~isfield(cnfg,'pval'), cnfg.pval =   0.05; end
%if ~isfield(cnfg,'pval_clus'), cnfg.pval_clus = '0.01'; end
if ~isfield(cnfg,'Nsurro'), cnfg.Nsurro = 0; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if isfield(cnfg,'outpath')
    if ~strcmp(cnfg.outpath(end),'\')
        cnfg.outpath = [cnfg.outpath '\'];
    end
end
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end
cnfg.Fs = ftdata.fsample;

%% COMPUTE GC

if iscell(cnfg.channel)
    Niter = length(cnfg.channel);
    channel = cnfg.channel;
else
    Niter = 1;
    channel{1} = cnfg.channel; %to have a cell format
end

for iter=1:Niter
    % Select data
    cfg = [];
    cfg.latency = cnfg.latency;
    cfg.channel = channel{iter};
    ftdata_aux      = ft_selectdata(cfg,ftdata);
    conc_data   = concatenate_fttrials(ftdata_aux);
    
    Fs = conc_data.fsample;
    if isfield(cnfg,'resample'), Fs = cnfg.resample; end
    data = conc_data.trial{1};
    if isfield(cnfg,'resample')
        for chi=1:size(data,1)
            data_aux(chi,:) = resample(data(chi,:),cnfg.resample,round(cnfg.Fs));
        end
        data=data_aux;
        clear data_aux
    end
    %data=zscore(data')';
    cfg = [];
    cfg.window   = round(cnfg.window*Fs);
    cfg.step     = round(cnfg.step*Fs);
    cfg.p        = cnfg.model_order;
    cfg.Nsurro   = cnfg.Nsurro; 
    cfg.Fs       = Fs;
    cfg.pval     = cnfg.pval;
    cfg.doplot   = cnfg.doplot;
    cfg.dosave   = cnfg.dosave;
    cfg.outpath  = cnfg.outpath;
    cfg.infosave = cnfg.infosave;
    cfg.label    = conc_data.label;
    if Niter>1
        cfg.infosave = [cfg.infosave '_case_' num2str(iter)];
    end
    cfg.mode     = 3; % Both temporal and freq. GC 
    
    data = diff(data')'; %Maybe this avoids some problems of stability
    [Favg{iter},favg{iter}]=computeGC_MVGC(data,cfg);
end

% GC = cell(1,size(cnfg.channel,1));
% 
% ch = cnfg.channel;
% f_env_v = (cnfg.f_env(1):cnfg.f_env(3):cnfg.f_env(2));
% BW    = cnfg.f_env(4);
% for chi = 1:size(ch,1)
%     data_x_aux = conc_data.trial{1}(ch(chi,1),:);
%     data_y_aux = conc_data.trial{1}(ch(chi,2),:);
% 
%     GC_xy = zeros(length(f_env_v),length(f_env_v));
%     GC_yx = zeros(length(f_env_v),length(f_env_v));
%     if cnfg.Nsurro>0
%         GC_xy_pval = zeros(length(f_env_v),length(f_env_v),cnfg.Nsurro);
%         GC_yx_pval = zeros(length(f_env_v),length(f_env_v),cnfg.Nsurro);
%     end
%     for x=1:length(f_env_v)
%         x_min = f_env_v(x)-BW/2;
%         x_max = f_env_v(x)+BW/2;
%         data_x = abs(hilbert(eegfilt(data_x_aux,cnfg.Fs,x_min,x_max)));
%         if isfield(cnfg,'resample')
%                 data_x = resample(data_x,cnfg.resample,cnfg.Fs); 
%         end
%         for y=1:length(f_env_v)
%             y_min = f_env_v(y)-BW/2;
%             y_max = f_env_v(y)+BW/2;
%             data_y = abs(hilbert(eegfilt(data_y_aux,cnfg.Fs,y_min,y_max)));
%             if isfield(cnfg,'resample')
%                 data_y = resample(data_y,cnfg.resample,cnfg.Fs); 
%             end
% 
%             %% GC HERE
%             [Favg,favg]=computeGC_MVGC([data_x' data_y']',cfg);
%             GC_xy(x,y) = Favg(2,1); %F(2,1) = GC(1->2)
%             GC_yx(x,y) = Favg(1,2);
%             %GCf_xy(x,y,:) = favg(2,1,:);
%             %GCf_yx(x,y,:) = favg(1,2,:);
%             %PSI_pval(x,y,:) = PSI_aux.pval;
%         end
%     end
%     GC{chi}.GC_xy = GC_xy;
%     GC{chi}.GC_yx = GC_yx;
%     %GC{chi}.GCf_xy = GCf_xy;
%     %GC{chi}.GCf_yx = GCf_yx;
%     GC{chi}.label_x = [conc_data.label{ch(chi,1)} '_env'];
%     GC{chi}.label_y = [conc_data.label{ch(chi,2)} '_env'];
%     GC{chi}.f_env_v = f_env_v;
%     % GC{chi}.pval --- I need to think about this
% end
% 
% 
% %% Plot and save results
% 
% if cnfg.dosave
%     if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
%     if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
% end
% 
% % One figure for each pair of channels
% if cnfg.doplot || cnfg.dosave
%     for iter=1:length(GC)
%         h=figure;
%         subplot(1,2,1)
%         imagesc([GC{iter}.f_env_v(1) GC{iter}.f_env_v(end)],...
%                 [GC{iter}.f_env_v(1) GC{iter}.f_env_v(end)],...
%                  GC{iter}.GC_xy);
%         title(['GC ' GC{iter}.label_x ' --> ' GC{iter}.label_y])
%         ylabel(['Freq. ' GC{iter}.label_x ' Hz'])
%         xlabel(['Freq. ' GC{iter}.label_y ' Hz'])
%         clim([0 max([GC{iter}.GC_xy(:)' GC{iter}.GC_yx(:)'])])
%         axis('xy')
%         subplot(1,2,2)
%         imagesc([GC{iter}.f_env_v(1) GC{iter}.f_env_v(end)],...
%                 [GC{iter}.f_env_v(1) GC{iter}.f_env_v(end)],...
%                  GC{iter}.GC_yx);
%         title(['GC ' GC{iter}.label_y ' --> ' GC{iter}.label_x])
%         ylabel(['Freq. ' GC{iter}.label_x ' Hz'])
%         xlabel(['Freq. ' GC{iter}.label_y ' Hz'])
%         clim([0 max([GC{iter}.GC_xy(:)' GC{iter}.GC_yx(:)'])])
%         axis('xy')
% 
%         pair = [GC{iter}.label_x '_' GC{iter}.label_y];
% 
%         if cnfg.dosave
%             savefig(h,[cnfg.outpath 'GC_' pair cnfg.infosave ])
%             saveas(h,[cnfg.outpath 'GC_' pair cnfg.infosave '.png'])
%         end
%         if ~cnfg.doplot, close(h), end
%     end
% end
% 
% % Save data
% if cnfg.dosave
%     save([cnfg.outpath 'GC' cnfg.infosave],'GC')
% end
% 
% end
% 
%  function newFigure1(h1,~)
% %% Function to act on subplot click
% %         Mouse click: Plots the selected subplot to a new figure
% %  Ctrl + Mouse click: Delete subplot
% 
%     switch get(gcf,'SelectionType')
%         case 'normal'
%             F = figure();
%             copyobj(h1.Children,gca(F));
%             % Copy the selected subplot title
%             tmp = get(h1,'title'); tmp = tmp.String;
%             % Set title to the new figure
%             title(gca(F), tmp);
%             
%             % Set axis to the new figure
%             tmp = get(h1,'XLim');
%             xlim(gca(F), tmp)
%             tmp = get(h1,'YLim');
%             ylim(gca(F), tmp)
%         case 'alt'
%             delete(h1);
%     end
%  end