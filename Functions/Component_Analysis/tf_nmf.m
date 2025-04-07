function [Wtf,H,W] = tf_nmf(cnfg,ftdata)

%% From continuos data (FieldTrip or matrix), compute Gabor, trials and NMF.
% It will do the trial average before the NMF.
% Time and Frequency are vectorize into a single dimension
% The second dimension is space (sensor)
% Therefore, it will compute:
%
%   X ~= WH
%
% - where X is the Gabor transform with time-freq in one dimension and
% sensors (space) in the second dimension
% - W are the tf distribution of each component
% - H are the spatial distribution of each component
%
% Syntax:
%    [Wtf,H,W] = tf_nmf(cfg,ftdata)
%
% Inputs:
%    cfg - Structure of parameters:
%
%       Ncomp      - Number of components 'k'. It can be a value or it can
%                    be estimated ('all','velicer')
%       stimdef    - [start end] Trial in seconds. Ex [-0.2 0.5]
%       time_trial - vector with each trial onset in seconds
%       M          - Time length for Gabor im samples (Def 128)
%       a          - Gabor resolution (Def M/16)
%       freqlim    - Frequencies of interest. Ex [5 100]
%       Nsurro     - Number of surrogates. Def=100
%       stats      - method for stats: Def 'std', 'cluster', 'pixel'
%       alpha      - pval for std stats. Def= 0.01
%
%       dosave   - logical. True/false save/not save the results
%       outpath  - string. Path to save the results if dosave=true
%       doplot   - logical. Plot the figure with the results.
%       infosave - string to include in the saved filed
%
%   data - Data in format field trip or matrix.
%
% %%%%%%%%%%% CASE 'data' is a matrix %%%%%%%%%%%%%%%%
% This code has been prepared to work with FieldTrip structures. However,
% it can be used with any matrix of data by defining some extra parameters:
%
%   data        - [Nch x Nsamples]
%   cfg.time    - vector of Nsamples with the timestamp of each sample
%   cfg.Fs      - Sampling rate
%   cfg.label   - [Optional] cell of strings with the name of each channel
%
% Outputs:
%    Wtf - W matrix reshaped to be 2D (time x freq)
%    H   - Spatial distribution
%    W   - TimeFreq distribution (vectorized)
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug. 2025; Last revision: 07-Apr-2025


%% Initialization

if ~isfield(cnfg,'Ncomp'), cnfg.Ncomp = 'velicer'; end
if ~isfield(cnfg,'M'), cnfg.M = 128; end
if ~isfield(cnfg,'a'), cnfg.a = cnfg.M/16; end
if ~isfield(cnfg,'Nsurro'), cnfg.Nsurro = 100; end
if ~isfield(cnfg,'stats'), cnfg.stats = 'std'; end
if ismember(cnfg.stats,{'cluster','pixel'}) && cnfg.Nsurro==0
    error(['You need to specify some surrogates to perform ' cnfg.stats ' stats'])
end
if ~isfield(cnfg,'alpha'), cnfg.alpha = 0.01; end
alpha = norminv(1-cnfg.alpha);

if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if ~isfield(cnfg,'outpath') && cnfg.dosave
    error('Outpath has not been specified to save the results'), end

%%% Case ftdata is a matrix
if ismatrix(ftdata)
    [Nch,Nsamples] = size(ftdata);
    disp(['Input data is a matrix with ' num2str(Nch) ' channels and '...
        num2str(Nsamples) ' samples.'])
    if ~isfield(cnfg,'time') 
        error('Parameter ''time'' is mandatory when working with matrix')
    end
    if ~isfield(cnfg,'Fs') 
        error('Parameter ''Fs'' is mandatory when working with matrix')
    end
    if ~isfield(cnfg,'label'), cnfg.label = []; end
    
    data = ftdata;
    clear ftdata
    ftdata.trial{1} = data;
    ftdata.fsample  = cnfg.Fs;
    ftdata.time{1}  = cnfg.time;
    ftdata.label    = cnfg.label;
end

%% COMPUTE GABOR

% Parameters
M = cnfg.M; %Time length
a = cnfg.a; % I divide the Fs by this 'a' value
time_trial = cnfg.time_trial;

%%%% GABOR TRANSFORM %%%%
disp('Computing GABOR Transform...')
y = ftdata.trial{1}';
Fs = ftdata.fsample;
g = gabwin({'tight', 'hann'}, a, M);
c = dgtreal(y, g, a, M);
ff = linspace(0,Fs/2,size(c,1));
pow_corrected = [];
for si=1:size(c,3)
    for fi=1:length(ff)
        pow_corrected(fi,:,si) = zscore(abs(c(fi,:,si)));
    end
end

% Limit the frequency dimension
if isfield(cnfg,'freqlim')
    [~, fini] = min(abs(ff - cnfg.freqlim(1)));  % index of initial freq
    [~, fend] = min(abs(ff - cnfg.freqlim(2)));  % index of initial freq    
    pow_corrected = pow_corrected(fini:fend,:,:);
end

% Create trials from continuous GABOR. Here I need time_trial
Ntrial = length(time_trial);
Fs_gabor = round(ftdata.fsample/a);
time_gabor = linspace(ftdata.time{1}(1),ftdata.time{1}(end),size(pow_corrected,2));
% Find closest index and value for each y
onset_samples = zeros(1,Ntrial);
for tri = 1:Ntrial
    [~, onset_samples(tri)] = min(abs(time_gabor - time_trial(tri)));  % index of minimum difference
end
pre_samples   = round(cnfg.stimdef(1)*Fs_gabor); % Samples before onset in samples
post_samples  = round(cnfg.stimdef(2)*Fs_gabor); % Samples after onset in samples
% Define trials
trl_ini = onset_samples + pre_samples;  % Start sample
trl_end = onset_samples + post_samples; % End sample
pow_trial_zscore = [];
for tri=1:Ntrial
    pow_trial_zscore(:,:,:,tri) = pow_corrected(:,trl_ini(tri):trl_end(tri),:);
end

pow_mean = mean(pow_trial_zscore,4);

%% NMF

%%% Vectorize
[T, F, S] = size(pow_mean);
pow_vec = reshape(pow_mean, [T*F, S]); 
V = pow_vec - min(pow_vec(:));

%%% Number of components
if isnumeric(cnfg.Ncomp)
    k = cnfg.Ncomp;
elseif strcmp(cnfg.Ncomp,'all')
    k = size(V,2);
elseif strcmp(cnfg.Ncomp,'velicer')
    [~, k] = velicer_map_oconnor(V);
end

% SVD-based initialization for stability
[U,D,S] = svd(V,"econ");
S = S';
U = U*D;
W0 = abs(U(:,1:k));
H0 = abs(S(1:k,:));

[W, H] = nnmf(V, k, 'w0', W0, 'h0', H0, 'algorithm', 'mult');

%%% Reshape to original dimensions
[T, F, ~] = size(pow_mean);
Wtf = reshape(W, [T, F, k]);   % Back to time × freq × sensors

%% STATS

if ismember(cnfg.stats,{'cluster','pixel'})
    % Values to keep later
    surro_cluster = zeros(1,cnfg.Nsurro);
    pixel_max = zeros(1,cnfg.Nsurro);
    pixel_min = zeros(1,cnfg.Nsurro);
    THclus_mask = zeros(T,F,k);
    THmax = [];
    
    disp('Computing surrogates on NMF...')
    k_surro=1;
    W0 = abs(U(:,1:k_surro));
    H0 = abs(S(1:k_surro,:));

    for si=1:cnfg.Nsurro
        % Define surrogate trials between the first and the last one at random time
        % points
        onset_ini = min(onset_samples);
        onset_end = max(onset_samples);
        onset_aux = randperm(onset_end-onset_ini);
        onset_surro = onset_aux(1:Ntrial)+onset_ini;
        trl_ini = onset_surro + pre_samples;  % Start sample
        trl_end = onset_surro + post_samples; % End sample
        pow_trial_surro = [];
        for tri=1:Ntrial
            pow_trial_surro(:,:,:,tri) = pow_corrected(:,trl_ini(tri):trl_end(tri),:);
        end
        
        pow_mean_surro = mean(pow_trial_surro,4);
        
        %%% Vectorize
        [Ts, Fs, Ss] = size(pow_mean_surro);
        pow_vec = reshape(pow_mean_surro, [Ts*Fs, Ss]);
        V = pow_vec - min(pow_vec(:));
        [W_surro, ~] = nnmf(V, k_surro, 'w0', W0, 'h0', H0, 'algorithm', 'mult');
        
        %%% For Pixel-based stats %%%
        pixel_max(si) = max(W_surro(:));
        pixel_min(si) = min(W_surro(:));
        
        %%% For cluster-stats %%%
        Wtf_surro = reshape(W_surro, [Ts, Fs, 1]);  
        % I put a thresholf of twice the std
        th_surro = mean(W_surro)+2*std(W_surro);
        Wtf_mask = Wtf_surro>=th_surro;
        % By default I select the maximum value as 'highest cluster'
        cluster_aux = max(W_surro(:));
        % Then, I do a loop over all clusters and select its value as the
        % sumatory of the values inside the cluster
        CC = bwconncomp(Wtf_mask);
        SS = regionprops(CC,'Area','PixelIdxList');
        for cli = 1:length(SS)
            cluster_aux = [cluster_aux sum(Wtf_surro(SS(cli).PixelIdxList))];
        end
        surro_cluster(si) = max(cluster_aux);
        
    end
      
    if strcmp(cnfg.stats,'pixel')
        % For each surrogate, I keep the pixel with highest and lowest weight
        % I define the THs as the 97.5% highest value in pixel_max and 97.5% lowest
        % value in pixel_min.
        % For now, only positive values (95%).   
        pixel_aux = sort(pixel_max);
        THmax = pixel_aux(round(cnfg.Nsurro*0.95));
        %pixel_aux = sort(pixel_min);
        %pixel_aux = pixel_aux(end:-1:1);
        %THmin = pixel_aux(round(cnfg.Nsurro*0.975)); 
        disp(['Threshold of significance pixel-based = ' num2str(THmax)])
    end
    
    if strcmp(cnfg.stats,'cluster')
        % For each surrogate, I keep the biggest cluster (sumatory of values inside)
        % I define the THs as the 97.5% highest value in pixel_max and 97.5% lowest
        % value in pixel_min.
        pixel_aux = sort(surro_cluster);
        THcluster = pixel_aux(round(cnfg.Nsurro*0.95));
        
        for ki=1:k
            % Cluster Original for each component
            Wtf_aux = squeeze(Wtf(:,:,ki));
            % I put a thresholf of twice the std
            th_orig = mean(Wtf_aux(:))+2*std(Wtf_aux(:));
            Wtf_mask = Wtf_aux>=th_orig;
            % Remove those clusters that are not significant
            CC = bwconncomp(Wtf_mask);
            SS = regionprops(CC,'Area','PixelIdxList');
            flag = [];
            for cli = 1:length(SS)
                if sum(Wtf_surro(SS(cli).PixelIdxList)) <= THcluster
                    flag = [flag cli];
                end
            end
            SS(flag)=[];
            % Create a mask of significance
            THclus_mask_aux = Wtf_aux*0;
            for cli = 1:length(SS)
                THclus_mask_aux(SS(cli).PixelIdxList) = 1;
            end
            THclus_mask(:,:,ki) = THclus_mask_aux;
        end
    end
end

%% PLOT results

if cnfg.dosave
    if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
    if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
end

if isfield(cnfg,'freqlim')
    ff = linspace(cnfg.freqlim(1),cnfg.freqlim(2),size(Wtf,1));
else
    ff = linspace(0,Fs/2,size(c,1));
end
tt = linspace(cnfg.stimdef(1),cnfg.stimdef(2),size(pow_mean,2));
label = ftdata.label;
    
if cnfg.doplot || cnfg.dosave
    
% Plot Time-Frequency: Wtf
h=figure;
for iter=1:k
    nrows = ceil(sqrt(k));
    ncols = ceil(k/nrows);
    G(iter) = subplot(nrows,ncols,iter); 
    G(iter).ButtonDownFcn = @newFigure1;
    hold on,
    Witer = squeeze(Wtf(:,:,iter));
    P=imagesc(tt,ff,Witer);
    
    %%% STATS %%%
    if strcmp(cnfg.stats,'pixel')
        mask_max = Witer>=THmax;
        % Check if image_processing_toolbox is installed
        if isempty(which('bwconncomp'))
            % Use my custom code
            warning('Using custom_bwboundaries, which has not been fully checked')
            B = custom_bwboundaries(mask_max);
        else
            B = bwboundaries(mask_max);
        end
        for i=1:length(B)
            plot(tt(B{i}(:,2)),ff(B{i}(:,1)),'k','LineWidth',1.5),
        end
        %mask_min = Witer<THmin;
        %B = bwboundaries(mask_min);
        %for i=1:length(B)
        %    plot(tt(B{i}(:,2)),ff(B{i}(:,1)),'k'),
        %end
    end
    if strcmp(cnfg.stats,'cluster')
        mask_max = squeeze(THclus_mask(:,:,iter));
        % Check if image_processing_toolbox is installed
        if isempty(which('bwconncomp'))
            % Use my custom code
            warning('Using custom_bwboundaries, which has not been fully checked')
            B = custom_bwboundaries(mask_max);
        else
            B = bwboundaries(mask_max);
        end
        for i=1:length(B)
            plot(tt(B{i}(:,2)),ff(B{i}(:,1)),'k','LineWidth',1.5),
        end
    end
    if strcmp(cnfg.stats,'std')
        Wtf_iter = squeeze(Wtf(:,:,iter));
        th_aux(1) = mean(Wtf_iter(:)+alpha*std(Wtf_iter(:)));
        th_aux(2) = mean(Wtf_iter(:)-alpha*std(Wtf_iter(:)));
        mask_max = Wtf_iter>=th_aux(1) | Wtf_iter<=th_aux(2);
        % Check if image_processing_toolbox is installed
        if isempty(which('bwconncomp'))
            % Use my custom code
            warning('Using custom_bwboundaries, which has not been fully checked')
            B = custom_bwboundaries(mask_max);
        else
            B = bwboundaries(mask_max);
        end
        for i=1:length(B)
            plot(tt(B{i}(:,2)),ff(B{i}(:,1)),'k','LineWidth',1.5),
        end
    end
    
    %colormap 'winter'
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    set(P,'UserData',iter);
    set(P, 'HitTest', 'off');
    axis xy
    axis([tt(1) tt(end) ff(1) ff(end)])
    %caxis([-max(abs(itc.tf_corrected(:))) max(abs(itc.tf_corrected(:)))])
    title(['TF Comp - ' num2str(iter)]); 
    colorbar
end
% Maximize the figure window to fill the screen
%set(h, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    if cnfg.dosave
        savefig(h,[cnfg.outpath 'NMF_Wtf' cnfg.infosave ])
        saveas(h,[cnfg.outpath 'NMF_Wtf' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h), end
end

if cnfg.dosave
    save([cnfg.outpath 'NMFval' cnfg.infosave],'Wtf','H','F','tt','ff','label')
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





