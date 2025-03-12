
function itc=ITPC_MEG_SEEG(cnfg,dataMEG,dataSEEG)

%% Computes inter-trial coherence between MEG and SEEG
%
%   TF = (TF - Baseline)/Baseline
%
% Time metrics: *not finished* 
%  - Inter-trial time-based correlation (ITCOR) *not finished*
%  - Inter-trial time-based partial correlation (pITCOR) *not finished*
%
% Time-frequency metrics:
%  - Inter-trial phase-based coherence (ITPC) *with surrogates*
%  - Inter-trial power-based coherence (ITLC) *with surrogates*
%  - Inter-site phase-based coherence (ISPC)  *with surrogates*
%  - Inter-site power-based coherence (ISPC)  *with surrogates*
%
% Syntax:  
%    itc=ITPC_MEG_SEEG(cfg,dataMEG,dataSEEG)
%
% Inputs:
%
%   cfg.channel     -
%      .channelMEG  -
%      .channelSEEG -
%      .metric      - {'TF', 'ITPC', 'ITLC', 'ISPC', 'ISLC'} 
%      .toi         - for example [0:0.01:1] - Only if trial data (no cont)
%      .foilim      - for example [1 45]
%      .trl         - [If cont data] output of ft_definetrial
%      .timetrl     - [If cont data] time of one trial
%      .trigger     - [Optional] for example [528 544] 
%      .stats       - ['fdr','surro']
%      .Nsurro      - [Optional] Default = 0
%      .FDR_q       - [Optional] Default = 0.05
%      .pval        - [Optional] Default = 0.01
%      .doplot      - true/false
%      .dosave      - true/false
%      .outpath     - string with path 
%      .infosave    - string with filename
%
% Outputs:
%
% See also:

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Sep. 2020; Last revision: 21-Aug-2024

% Change log:
% 21/08/2024: Include stats with FDR
% 18/12/2023: Updated to accept outpout from timelockanalysis
% 13/08/2021: Included the option to plot the time-freq map

% TO DO:
%   - Include option to compute tf on cont data and then create trials
%   - ITCOR and pITCOR
%   - Modify options baseline
%   - Separate or combine triggers (now it combines them)
% 

%% Default values

if ~isfield(cnfg,'channel')
    cnfg.channel  = 'all'; end
if ~isfield(cnfg,'channelMEG')
    cnfg.channelMEG  = 'all'; end
if ~isfield(cnfg,'channelSEEG')
    cnfg.channelSEEG = 'all'; end
if ~isfield(cnfg,'toi') && ~isfield(cnfg,'timetrial')
    error('The field toi or timetrial are mandatory'); end
if isfield(cnfg,'timetrial')
    cnfg.toi = cnfg.timetrial; end
if ~isfield(cnfg,'foilim')
    error('The field foilim is mandatory'); end
if ~isfield(cnfg,'baseline')
    error('The field baseline is mandatory'); end
if ~isfield(cnfg,'Nsurro')
    cnfg.Nsurro=0; end
if ~isfield(cnfg,'trigger')
    if isfield(dataMEG,'trialinfo')
        cnfg.trigger=unique(dataMEG.trialinfo); end
end
if ~isfield(cnfg,'stats'), cnfg.stats='fdr'; end
if ~isfield(cnfg,'FDR_q')
    cnfg.FDR_q=0.05; end
if ~isfield(cnfg,'pval')
    cnfg.pval=0.01; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if isfield(cnfg,'outpath')
    if ~strcmp(cnfg.outpath(end),'\')
        cnfg.outpath = [cnfg.outpath '\'];
    end
end
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end

if ~isfield(cnfg,'metric')
    cnfg.metric = {'ITPC', 'ITLC', 'ISPC', 'ISPC'}; end
if ischar(cnfg.metric)
    metric = cnfg.metric;
    cnfg.metric = [];
    cnfg.metric{1} = metric; 
end
m = ismember({'TF', 'ITPC', 'ITLC', 'ISPC', 'ISLC'}, cnfg.metric);
if sum(m)==0
    error('No metrics have been selected.'); end

if nargin == 2 %Not MEG-SEEG mix
    if isfield(cnfg,'trigger')
    trials = find(ismember(dataMEG.trialinfo,cnfg.trigger));
    cfg = [];
    cfg.channel = cnfg.channel;
    if ~isfield(cnfg,'trl')
        cfg.trials = trials; end
    dataMEG = ft_selectdata(cfg,dataMEG);
    end
    dataMIX = dataMEG;
    
else %Create MEG-SEEG struct  
    trials = find(ismember(dataMEG.trialinfo,cnfg.trigger));
    cfg = [];
    cfg.channel = cnfg.channelMEG;
    cfg.trials = trials;
    dataMEG=ft_selectdata(cfg,dataMEG);
    cfg = [];
    cfg.channel = cnfg.channelSEEG;
    cfg.trials = trials;
    dataSEEG=ft_selectdata(cfg,dataSEEG);
    
    dataMIX.label = [dataMEG.label' dataSEEG.label']';
    dataMIX.time = dataMEG.time;
    dataMIX.fsample = dataMEG.fsample;
    for tr=1:length(dataMEG.trial)
        dataMIX.trial{1,tr} = [dataMEG.trial{tr}' dataSEEG.trial{tr}']';
    end
end

if ~strcmp(cnfg.stats,'surro')
    cnfg.Nsurro=0;
end

%% Compute original ITC

disp('Computing wavelets...')
cfg        = [];
cfg.method = 'wavelet';
cfg.output = 'fourier';
cfg.foilim = cnfg.foilim;
if isfield(cnfg,'trl') %cont data
    cfg.toi    = dataMIX.time{1};
else
    cfg.toi    = cnfg.toi; %trial data
end
cfg.pad='nextpow2';
freq       = ft_freqanalysis(cfg, dataMIX);
freq.fourierspctrm(isnan(freq.fourierspctrm))=0;
if isfield(cnfg,'trl') %cont data
    freq = create_trials(freq,cnfg); % I need to finish this
end

% if ismember('TF', cnfg.metric)
%     cfg.output = 'pow';
%     freq_power = ft_freqanalysis(cfg, dataMIX);
% end
freq.powspctrm = squeeze(mean(abs(freq.fourierspctrm.^2),1));
if length(freq.label)==1 %Only one channel
    pow_aux = freq.powspctrm;
    freq.powspctrm = [];
    freq.powspctrm(1,:,:,:) = pow_aux;
end

% make a new FieldTrip-style data structure containing the ITC
% copy the descriptive fields over from the frequency decomposition

itc           = [];
itc.label     = freq.label;
itc.freq      = freq.freq;
itc.time      = freq.time;
itc.dimord    = 'chan_freq_time';
itc.cfg       = cnfg;

F = freq.fourierspctrm;   % copy the Fourier spectrum
N = size(F,1);           % number of trials

% compute inter-trial phase coherence (itpc)
if ismember('ITPC', cnfg.metric)
    disp('Computing inter-trial phase coherence...')
    itc.itpc=compute_ITPC(F);
    
%     %%% Surrogates
%     ntp = size(F,4);
%     %nfp = size(F,3);
%     %nch = size(F,2);
%     if cnfg.Nsurro>0
%     if size(F,1)>(ntp+4)
%         warning('More trials than time-points. Surrogates cannot be computed')
%     else
%     parfor s=1:cnfg.Nsurro
%         %Surro F
%         p=randperm(ntp-4); p=p+2;
%         F_surro = F*0;
%             
%         %This is wrong, but there must be a way to create idx. Mych
%         %faster
% %         for ntr=1:N
% %             Fidx_aux(ntr,:) = [p(ntr)+1:ntp 1:p(ntr)];
% %         end
% %         Fidx = repmat(Fidx_aux,1,1,nch,nfp);
% %         Fidx = permute(Fidx,[1 3 4 2]);
% %         F_surro = F(Fidx);
%             
%         for ntr=1:N
%             F_surro(ntr,:,:,:) = [squeeze(F(ntr,:,:,p(ntr)+1:end)) squeeze(F(ntr,:,:,1:p(ntr)))];
%         end
%         itpc_surro(:,:,:,s)=compute_ITPC(F_surro);
%     end
%     itc.itpc_surro = itpc_surro;
%     end
%     end
end

% compute inter-trial linear coherence (itlc)
if ismember('ITLC', cnfg.metric)
    disp('Computing inter-trial linear coherence...')
    itc.itlc      = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
    itc.itlc      = abs(itc.itlc);     % take the absolute value, i.e. ignore phase
    itc.itlc      = squeeze(itc.itlc); % remove the first singleton dimension
end

% compute inter-site phase coherence (ispc)
if ismember('ISPC', cnfg.metric)
    disp('Computing inter-site phase coherence...')
    Nch = size(F,2);
    c=0;
    for ni=1:Nch-1
        for nj=ni+1:Nch
            c=c+1;
            itc.label_pair{c,1} = [itc.label{ni} ' vs ' itc.label{nj}];
            itc.ispc.ispc(c,:,:) = squeeze(abs(mean(exp(1i*(angle(F(:,ni,:,:))-angle(F(:,nj,:,:)))))));
        end
    end
end

% "Time frequency inter-trial coherence"
% compute inter-site linear coherence (islc)
if ismember('ISLC', cnfg.metric)
    disp('Computing inter-site linear coherence...')
    Nf = size(F,3);
    Nt = size(F,4);
    Vt = (2*Nt*Nt+1) : 2*Nt+1 : 4*Nt*Nt;
    c=0;
    for ni=1:Nch-1
        for nj=ni+1:Nch
            c=c+1;
            for f=1:Nf 
                COR = corrcoef([abs(squeeze(F(:,ni,f,:))) abs(squeeze(F(:,nj,f,:)))]);
                itc.islc.islc(c,f,:) = COR(Vt);
            end
        end
    end
end

%% Statistical significance

% For TF, we test the difference between ERP and baseline in decibels.
% We can correct these values usins FDR or Surrogate analysis
% With surrogates, we test whether this difference is significative.

if ismember({'TF'}, cnfg.metric)
    t_bl(1) = find(freq.time>=cnfg.baseline(1),1);
    t_bl(2) = find(freq.time>=cnfg.baseline(2),1);
    %t0 = find(freq.time>=0,1);
    % Include also the baseline in the stats
    t0 = 1;
    
    %Ratio in dB between TF values versus baseline
    MU_power = mean(freq.powspctrm(:,:,t_bl(1):t_bl(2)),3);
    for t=1:size(freq.powspctrm,3)
        %freq.power_corrected(:,:,t) = 10*log10(freq.powspctrm(:,:,t)./MU_power);
        % Updated to have the relative ratio with respect to baseline
        freq.power_corrected(:,:,t) = (freq.powspctrm(:,:,t)-MU_power)./MU_power;
    end
    % Correct for non-measured values
    freq.power_corrected(freq.power_corrected==-1) = NaN;
    freq.powspctrm(freq.powspctrm==0) = NaN;
    freq.fourierspctrm(freq.fourierspctrm==0) = NaN;
    itc.powspctrm    = freq.powspctrm; 
    itc.tf_corrected = freq.power_corrected;

    %%% Stats using FDR
    % Compute pval for each point across trials.
    % Test if TF is different from baseline.
    [Ntr,Nch,Nf,Nt] = size(F); 
    pval_tf = ones(Nch,Nf,Nt);
    mask_tf = zeros(Nch,Nf,Nt);
    for chi=1:Nch
        baseline_aux = squeeze(mean(abs(freq.fourierspctrm(:,chi,:,t_bl(1):t_bl(2)).^2),4));
        for ti=t0:Nt
            power_aux = squeeze(abs(freq.fourierspctrm(:,chi,:,ti).^2));
            [~,pval_tf(chi,:,ti)] = ttest(baseline_aux,power_aux);
            %rel_ratio = (power_aux-baseline_aux)./baseline_aux;
            %[~,pval_tf(chi,:,ti)] = ttest(rel_ratio);
        end
        [pID,pN] = fdr(pval_tf(chi,:,t0:end),cnfg.FDR_q); %False Discovery Rate for each channel
        mask_tf(chi,:,:)  = pval_tf(chi,:,:)<pID;
        itc.fdr_tf(chi,:) = [pID,pN];
    end
    itc.pval_tf = pval_tf;
     
    %%% Stats using surrogates
    if cnfg.Nsurro>0 && strcmp(cnfg.stats,'surro')
        freq.power = abs(freq.fourierspctrm.^2);
        ntp = size(freq.power,4);
    if size(freq.power,1)>(ntp+4)
        warning('More trials than time-points. Surrogates cannot be computed')
    else
    disp('Computing surrogates for time-frequency analysis...')
    freq.power = permute(freq.power,[2 4 3 1]);
    for s=1:cnfg.Nsurro
        %Surro F
        p=randperm(ntp-4); p=p+2;
        power_surro = freq.power*0;
        if length(freq.label)>1 %More than one channel        
            for ntr=1:N
                power_surro(:,:,:,ntr) = [squeeze(freq.power(:,p(ntr)+1:end,:,ntr)) squeeze(freq.power(:,1:p(ntr),:,ntr))];
            end
            power_surro = permute(power_surro,[4 1 3 2]);
            power_surro_mean = squeeze(mean(power_surro,1));
        else
            for ntr=1:N
                power_surro(:,:,:,ntr) = vertcat(squeeze(freq.power(:,p(ntr)+1:end,:,ntr)),squeeze(freq.power(:,1:p(ntr),:,ntr)));
            end            
            power_surro = permute(power_surro,[4 1 3 2]);
            power_surro_mean = squeeze(mean(power_surro,1));
            pow_aux = power_surro_mean;
            power_surro_mean = [];
            power_surro_mean(1,:,:,:) = pow_aux;
        end

        MU_surro = mean(power_surro_mean(:,:,t_bl(1):t_bl(2)),3);
        for t=1:size(power_surro_mean,3)
            power_surro_corrected(:,:,t,s) = 10*log10(power_surro_mean(:,:,t)./MU_surro);
        end
    end
    freq.power = permute(freq.power,[4 1 3 2]);
    %%% Pixel-based stats
    for ch=1:size(power_surro_corrected,1)
        for s=1:cnfg.Nsurro    
            power_aux = power_surro_corrected(ch,:,:,s);
            power_pixel(s) = max(power_aux(:));
            power_pixel_min(s) = min(power_aux(:));
        end
        itc.tf_th(ch) = quantile(power_pixel,1-cnfg.pval);
        itc.tf_th_inf(ch) = quantile(power_pixel_min,cnfg.pval);
    end
    end 
    end
end

% For ITPC, first we correct computing the p-val with respect to the
% baseline. Then, if surrogates were computed, we can correct for multiple
% comparisons using cluster-based stats.

if isfield(itc,'itpc')
    t_bl(1) = find(freq.time>=cnfg.baseline(1),1);
    t_bl(2) = find(freq.time>=cnfg.baseline(2),1);
    
    %Difference between ITPC values versus baseline
    MU_itpc = mean(itc.itpc(:,:,t_bl(1):t_bl(2)),3);
    for t=1:size(itc.itpc,3)
        itpc_diff(:,:,t) = (itc.itpc(:,:,t) - MU_itpc);
    end
    
    % Pvalues using Rayleigh’s Z (Cohen 2014)
    itc.itpc_pval=exp(-N*itpc_diff.^2);
    mask_itpc = itc.itpc_pval<=cnfg.pval & itpc_diff>0;
    itc.itpc_sig = itc.itpc.*mask_itpc;
    
    %%% TO DO %%%
%     % Correct for multple comparisons
%     if isfield(itc,'itpc_surro')
%     end
end

% We follow a surrogate approach or inter-site coherence, randomly
% alternating the trials. Then, we compute the p-value at each pixel and
% correct for multiple comparisons using FDR.

m = ismember({'ISPC', 'ISLC'}, cnfg.metric);
if sum(m) > 0
Ns = cnfg.Nsurro;
if Ns>0
    q = cnfg.FDR_q;
    disp(['Computing significance using ' num2str(Ns) ' surrogates...'])
    
    Nch = size(F,2);
    barra = waitbar(0,['Computing surrogates: 0/' num2str(Ns)]);
    for s=1:Ns
        %Randomize trials
        F_surro = F;
        for n=1:Nch
            F_surro(:,n,:,:) = F(randperm(size(F,1)),n,:,:);
        end
        
        
        % compute inter-site phase coherence (ispc)
        if ismember('ISPC', cnfg.metric)
        c=0;
        for ni=1:Nch-1
            for nj=ni+1:Nch
                c=c+1;
                itc.ispc.ispc_surro(c,:,:,s) = squeeze(abs(mean(exp(1i*(angle(F_surro(:,ni,:,:))-angle(F_surro(:,nj,:,:)))))));
            end
        end
        end
        
        % compute inter-site linear coherence (islc)
        if ismember('ISLC', cnfg.metric)
        c=0;
        for ni=1:Nch-1
            for nj=ni+1:Nch
                c=c+1;
                for f=1:Nf
                    COR = corrcoef([abs(squeeze(F_surro(:,ni,f,:))) abs(squeeze(F_surro(:,nj,f,:)))]);
                    itc.islc.islc_surro(c,f,:,s) = COR(Vt);
                end
            end
        end
        end
        
        barra=waitbar(s/Ns,barra,['Computing surrogates: ' num2str(s) '/' num2str(Ns)]);
    end
    close(barra)
    
    if ismember('ISPC', cnfg.metric)
        ispc_mean = mean(itc.ispc.ispc_surro,4);
        ispc_std  = permute(std(permute(itc.ispc.ispc_surro,[4 1 2 3])),[2 3 4 1]);
        itc.ispc.ispc_pval = normcdf(itc.ispc.ispc,ispc_mean,ispc_std);
        itc.ispc.ispc_pval = 1-itc.ispc.ispc_pval;
        for iter = 1:size(itc.ispc.ispc_pval)
            [pID,pN] = fdr(itc.ispc.ispc_pval(iter,:,:),q);
            if pID>=0, itc.ispc.pID(iter)=pID; else, itc.ispc.pID(iter)=0; end
            if pN>=0, itc.ispc.pN(iter)=pN; else, itc.ispc.pN(iter)=0; end
        end
        
        itc.ispc.ispc_surro = []; %Delete it to save space. Save mean and std instead
        itc.ispc.ispc_mean_surro = ispc_mean;
        itc.ispc.ispc_std_surro = ispc_std;
    end
    
    if ismember('ISLC', cnfg.metric)
        islc_mean = mean(itc.islc.islc_surro,4);
        islc_std  = permute(std(permute(itc.islc.islc_surro,[4 1 2 3])),[2 3 4 1]);
        itc.islc.islc_pval = normcdf(itc.islc.islc,islc_mean,islc_std);
        itc.islc.islc_pval = 1-itc.islc.islc_pval;
        for iter = 1:size(itc.islc.islc_pval)
            [pID,pN]=fdr(itc.islc.islc_pval(iter,:,:),q);
            if pID>=0, itc.islc.pID(iter)=pID; else, itc.islc.pID(iter)=0; end
            if pN>=0, itc.islc.pN(iter)=pN; else, itc.islc.pN(iter)=0; end
        end
        
        itc.islc.islc_surro = []; %Delete it to save space. Save mean and std instead
        itc.islc.islc_mean_surro = islc_mean;
        itc.islc.islc_std_surro = islc_std;
    end
end
end

%% PLOT results

if cnfg.dosave
    if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
    if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
end

if cnfg.doplot || cnfg.dosave
    
% Plot Time-Frequency
if ismember('TF', cnfg.metric)
h_tf=figure;
for iter=1:size(itc.tf_corrected,1)
    nrows = ceil(sqrt(length(freq.label)));
    ncols = ceil(length(freq.label)/nrows);
    G(iter) = subplot(nrows,ncols,iter); 
    G(iter).ButtonDownFcn = @newFigure1;
    
    hold on
    tf_fig = squeeze(itc.tf_corrected(iter,:,:));
    if isfield(itc,'pval_tf')
        tf_fig_pval = squeeze(itc.pval_tf(iter,:,:));
        %correct NaN values and zero values
        %tf_fig_pval(isnan(tf_fig_pval))=1;
        %tf_fig_pval(tf_fig_pval==0)=1;

        % Option 1: Apply a mask
        %tf_fig(tf_fig_pval>itc.fdr_tf(iter,1)) = 0; %Select here the FDR. 1=pID; 2=pN
        
        % Option 2: Plot all the TF and include a circle where it is
        % signficant
        mask_sig = tf_fig_pval<itc.fdr_tf(iter,1);
        B = bwboundaries(mask_sig);
        %tf(:,1:7,:)=0;
        %figure, hold on,
        %imagesc(itc.time, itc.freq, squeeze(mean(tf_osc,1))-squeeze(mean(tf_pt,1)));
    end
    if isfield(itc,'tf_th')
        tf_fig(tf_fig<itc.tf_th(iter) & tf_fig>itc.tf_th_inf(iter)) = 0;
    end
    P=imagesc(itc.time, itc.freq, tf_fig);
    % Option 2: Plot circle of significance
    for i=1:length(B)
        plot(itc.time(B{i}(:,2)),itc.freq(B{i}(:,1)),'k'),
    end
    
    %colormap 'winter'
    set(P,'UserData',iter);
    set(P, 'HitTest', 'off');
    axis xy
    axis([cnfg.toi(1) cnfg.toi(end) cnfg.foilim(1) cnfg.foilim(end)])
    caxis([-max(abs(itc.tf_corrected(:))) max(abs(itc.tf_corrected(:)))])
    title(['TF - ' itc.label{iter}]); 
    hold off
    
%     hold on
%     cfg              = [];
%     if isfield(cnfg,'baseline')
%         cfg.baseline     = cnfg.baseline; end
%     cfg.channel      = iter;
%     cfg.baselinetype = 'zscore';
%     ft_singleplotTFR(cfg, freq); 
%     %colormap 'winter'
%     title(['TF - ' itc.label{iter}]);  
%     hold off
end
% Maximize the figure window to fill the screen
set(h_tf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    if cnfg.dosave
        savefig(h_tf,[cnfg.outpath 'TF' cnfg.infosave ])
        saveas(h_tf,[cnfg.outpath 'TF' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h_tf), end
end


% Plot ITPC
if ismember('ITPC', cnfg.metric)
h_itpc=figure;
for iter=1:size(itc.itpc,1)
    nrows = ceil(sqrt(size(itc.itpc,1)));
    ncols = ceil(size(itc.itpc,1)/nrows);
    G(iter) = subplot(nrows,ncols,iter); 
    G(iter).ButtonDownFcn = @newFigure1;
    %%%%%%%%%%%%%%%%%%%%%%%%
    hold on
    if isfield(itc,'itpc_sig')
        % Option 1: Plot mask
        %P=imagesc(itc.time, itc.freq, squeeze(itc.itpc_sig(iter,:,:)));
        
        % Option 2: plot all data and include significance
        P=imagesc(itc.time, itc.freq, squeeze(itc.itpc(iter,:,:)));
        mask_sig = squeeze(itc.itpc_sig(iter,:,:))>0;
        B = bwboundaries(mask_sig);
        for i=1:length(B)
            plot(itc.time(B{i}(:,2)),itc.freq(B{i}(:,1)),'k'),
        end
    else
        P=imagesc(itc.time, itc.freq, squeeze(itc.itpc(iter,:,:)));
    end
    %colormap 'winter'
    set(P,'UserData',iter);
    set(P, 'HitTest', 'off');
    axis xy
    axis([cnfg.toi(1) cnfg.toi(end) cnfg.foilim(1) cnfg.foilim(end)])
    caxis([0 max(itc.itpc(:))])
    title(['ITPC - ' itc.label{iter}]);  
    hold off
end
% Maximize the figure window to fill the screen
set(h_itpc, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    if cnfg.dosave
        savefig(h_itpc,[cnfg.outpath 'ITPC' cnfg.infosave ])
        saveas(h_itpc,[cnfg.outpath 'ITPC' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h_itpc), end
end

% Plot ITLC
if ismember('ITLC', cnfg.metric)
h_itlc=figure;
for iter=1:size(itc.itlc,1)
    nrows = ceil(sqrt(size(itc.itlc,1)));
    ncols = ceil(size(itc.itlc,1)/nrows);
    G(iter) = subplot(nrows,ncols,iter); 
    G(iter).ButtonDownFcn = {@newFigure1};
    
    hold on,
    P=imagesc(itc.time, itc.freq, squeeze(itc.itlc(iter,:,:)));
    axis xy
    axis([cnfg.toi(1) cnfg.toi(end) cnfg.foilim(1) cnfg.foilim(end)])
    caxis([0 max(itc.itlc(:))])
    set(P,'UserData',iter);
    set(P, 'HitTest', 'off');
    title(['ITLC - ' itc.label{iter}]);   
    hold off
end
    if cnfg.dosave
        savefig(h_itlc,[cnfg.outpath 'ITLC' cnfg.infosave ])
        saveas(h_itlc,[cnfg.outpath 'ITLC' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h_itlc), end
end

% Plot ISPC
if ismember('ISPC', cnfg.metric)
h_ispc=figure;
for iter=1:size(itc.ispc.ispc,1)
    nrows = ceil(sqrt(size(itc.ispc.ispc,1)));
    ncols = ceil(size(itc.ispc.ispc,1)/nrows);
    G(iter) = subplot(nrows,ncols,iter); 
    G(iter).ButtonDownFcn = @newFigure1;
    
    hold on,
    P=imagesc(itc.time, itc.freq, squeeze(itc.ispc.ispc(iter,:,:)));
    axis xy
    axis([cnfg.toi(1) cnfg.toi(end) cnfg.foilim(1) cnfg.foilim(end)])
    set(P,'UserData',iter);
    set(P, 'HitTest', 'off');
    title(['ISPC - ' itc.label_pair{iter}]);  
end
    if cnfg.dosave
        savefig(h_ispc,[cnfg.outpath 'ISPC' cnfg.infosave ])
        saveas(h_ispc,[cnfg.outpath 'ISPC' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h_ispc), end
end

% Plot significant ISPC
if ismember('ISPC', cnfg.metric)
h_ispc_pval=figure;
for iter=1:size(itc.ispc.ispc,1)
    nrows = ceil(sqrt(size(itc.ispc.ispc,1)));
    ncols = ceil(size(itc.ispc.ispc,1)/nrows);
    G(iter) = subplot(nrows,ncols,iter); 
    G(iter).ButtonDownFcn = @newFigure1;
    
    hold on,
    ispc=squeeze(itc.ispc.ispc(iter,:,:));
    ispc_pval=squeeze(itc.ispc.ispc_pval(iter,:,:));
    ispc(ispc_pval>itc.ispc.pN(iter))=0;
    P=imagesc(itc.time, itc.freq, ispc);
    axis xy
    axis([cnfg.toi(1) cnfg.toi(end) cnfg.foilim(1) cnfg.foilim(end)])
    set(P,'UserData',iter);
    set(P, 'HitTest', 'off');
    title(['ISPC - ' itc.label_pair{iter}]);  
end
    if cnfg.dosave
        savefig(h_ispc_pval,[cnfg.outpath 'ISPCpval' cnfg.infosave ])
        saveas(h_ispc_pval,[cnfg.outpath 'ISPCpval' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h_ispc_pval), end
end

% Plot ISLC
if ismember('ISLC', cnfg.metric)
h_islc=figure;
for iter=1:size(itc.islc.islc,1)
    nrows = ceil(sqrt(size(itc.islc.islc,1)));
    ncols = ceil(size(itc.islc.islc,1)/nrows);
    G(iter) = subplot(nrows,ncols,iter); 
    G(iter).ButtonDownFcn = @newFigure1;
    
    hold on,
    P=imagesc(itc.time, itc.freq, squeeze(itc.islc.islc(iter,:,:)));
    axis xy
    axis([cnfg.toi(1) cnfg.toi(end) cnfg.foilim(1) cnfg.foilim(end)])
    set(P,'UserData',iter);
    set(P, 'HitTest', 'off');
    title(['ISLC - ' itc.label_pair{iter}]);   
    hold off
end
    if cnfg.dosave
        savefig(h_islc,[cnfg.outpath 'ISLC' cnfg.infosave ])
        saveas(h_islc,[cnfg.outpath 'ISLC' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h_islc), end
end

% Plot significant ISLC
if ismember('ISLC', cnfg.metric)
h_islc_pval=figure;
for iter=1:size(itc.islc.islc,1)
    nrows = ceil(sqrt(size(itc.islc.islc,1)));
    ncols = ceil(size(itc.islc.islc,1)/nrows);
    G(iter) = subplot(nrows,ncols,iter); 
    G(iter).ButtonDownFcn = @newFigure1;
    
    hold on,
    islc=squeeze(itc.islc.islc(iter,:,:));
    islc_pval=squeeze(itc.islc.islc_pval(iter,:,:));
    islc(islc_pval>itc.islc.pN(iter))=0;
    P=imagesc(itc.time, itc.freq, islc);
    axis xy
    axis([cnfg.toi(1) cnfg.toi(end) cnfg.foilim(1) cnfg.foilim(end)])
    set(P,'UserData',iter);
    set(P, 'HitTest', 'off');
    title(['ISLC - ' itc.label_pair{iter}]);  
end
    if cnfg.dosave
        savefig(h_islc_pval,[cnfg.outpath 'ISLCpval' cnfg.infosave ])
        saveas(h_islc_pval,[cnfg.outpath 'ISLCpval' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h_islc_pval), end
end
end

if cnfg.dosave
    %pairnodes = [ftdata.label{chi} '_' ftdata.label{chj+length(cnfg.SEEGch)}];
    %pairnodes = [wcoh.labelcmb{ch,2} '_' wcoh.labelcmb{ch,1}];
    %savefig(h,[cnfg.outpath 'WCOH_connectivity_' cnfg.infosave '_' pairnodes])
    %saveas(h,[cnfg.outpath 'WCOH_connectivity_' cnfg.infosave '_' pairnodes '.png'])
    save([cnfg.outpath 'ITCval' cnfg.infosave],'itc')
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

% %% compute inter-trial phase coherence (itpc)
% function itpc=compute_ITPC(F)
%     N = size(F,1); 
%     itpc = F./abs(F);         % divide by amplitude
%     itpc = sum(itpc,1);   % sum angles
%     itpc = abs(itpc)/N;   % take the absolute value and normalize
%     itpc = squeeze(itpc); % remove the first singleton dimension
%     
%     if size(F,2)==1 %Only one channel
%         itpc_aux = itpc;
%         itpc = [];
%         itpc(1,:,:) = itpc_aux;
%     end
% end

%% Create trials
function freq = create_trials(freq,cnfg)

end
