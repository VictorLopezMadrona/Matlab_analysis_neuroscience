function pow_mean=tf_gabor(cnfg,ftdata)

%% From continuos data (FieldTrip or matrix), compute Gabor and TF response.
%
% This code is much faster than ITPC_MEG_SEEG 
%
% Syntax:
%    tf_gabor(cnfg,ftdata)
%
% Inputs:
%    cfg - Structure of parameters:
%
%       stimdef    - [start end] Trial in seconds. Ex [-0.2 0.5]
%       time_trial - vector with each trial onset in seconds
%       M          - Time length for Gabor im samples (Def 128)
%       a          - Gabor resolution (Def M/16)
%       freqlim    - Frequencies of interest. Ex [5 100]
%       stats      - method for stats: Def 'std'
%       alpha      - pval for std stats. Def= 0.01
%       baseline   - [tmin tmax]; 'all' (def); 'none'
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
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug. 2025; Last revision: 14-Apr-2025


%% Initialization

if ~isfield(cnfg,'M'), cnfg.M = 128; end
if ~isfield(cnfg,'a'), cnfg.a = cnfg.M/16; end
if ~isfield(cnfg,'stats'), cnfg.stats = 'std'; end
if ~isfield(cnfg,'alpha'), cnfg.alpha = 0.01; end
alpha = norminv(1-cnfg.alpha);
if ~isfield(cnfg,'baseline'), cnfg.baseline = 'all'; end
if ~isfield(cnfg,'freqlim'), cnfg.freqlim = [0 ftdata.fsample/2]; end

if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if ~isfield(cnfg,'outpath') && cnfg.dosave
    error('Outpath has not been specified to save the results'), end

%%% Case ftdata is a matrix
if ~isstruct(ftdata)
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

if size(ftdata.trial{1},1) > 50
    warning('Too many signals. It may cause problems of memory') 
    warning('Try gabor_iter.m instead')
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

%%% Correct power - NO BASELINE SELECTION FOR NOW   
% Case zscore each frequency
if strcmp(cnfg.baseline,'all')
%     for fi=1:length(ff)
%         c_aux = abs(c(fi,:,:));
%         for si=1:size(c,3)
%             % Manual zscore on each channel using the info from all channels
%             pow_corrected(fi,:,si) = (abs(c(fi,:,si))-mean(c_aux(:))) / std(c_aux(:)) ;
%         end
%     end
    pow_corrected = [];
    for si=1:size(c,3)
       for fi=1:length(ff)
           pow_corrected(fi,:,si) = zscore(abs(c(fi,:,si)));
       end
    end
    
%No correction    
else
    pow_corrected = abs(c);
end

% Limit the frequency dimension
%if isfield(cnfg,'freqlim')
%    [~, fini] = min(abs(ff - cnfg.freqlim(1)));  % index of initial freq
%    [~, fend] = min(abs(ff - cnfg.freqlim(2)));  % index of initial freq    
%    pow_corrected = pow_corrected(fini:fend,:,:);
%end

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

%% PLOT results

if size(pow_mean,3)>50 
    warning('Too many signals. No plot not saving will be done')
    cnfg.dosave = false;
    cnfg.doplot = false;
end
    
if cnfg.dosave
    if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
    if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
end

%if isfield(cnfg,'freqlim')
%    ff = linspace(cnfg.freqlim(1),cnfg.freqlim(2),size(c,1));
%else
%    ff = linspace(0,Fs/2,size(c,1));
%end
tt = linspace(cnfg.stimdef(1),cnfg.stimdef(2),size(pow_mean,2));
label = ftdata.label;
    
if cnfg.doplot || cnfg.dosave
    
% Plot Time-Frequency
h=figure;
Nch = size(pow_mean,3);
for iter=1:Nch
    nrows = ceil(sqrt(Nch));
    ncols = ceil(Nch/nrows);
    G(iter) = subplot(nrows,ncols,iter); 
    G(iter).ButtonDownFcn = @newFigure1;
    hold on,
    Pw_iter = squeeze(pow_mean(:,:,iter));
    P=imagesc(tt,ff,Pw_iter);
    
    %%% STATS %%%
    if strcmp(cnfg.stats,'std')
        Wtf_iter = squeeze(pow_mean(:,:,iter));
        th_aux(1) = mean(Wtf_iter(:))+alpha*std(Wtf_iter(:));
        th_aux(2) = mean(Wtf_iter(:))-alpha*std(Wtf_iter(:));
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
    axis([tt(1) tt(end) cnfg.freqlim])
    title(['TF Comp - ' num2str(iter)]); 
    colorbar
end
% Maximize the figure window to fill the screen
set(h, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
if cnfg.dosave
    savefig(h,[cnfg.outpath 'TF_Gabor' cnfg.infosave ])
    saveas(h,[cnfg.outpath 'TF_Gabor' cnfg.infosave '.png'])
end
if ~cnfg.doplot, close(h), end
end

%if cnfg.dosave
%    save([cnfg.outpath 'NMFval' cnfg.infosave],'Wtf','H','F','tt','ff','label')
%end

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

