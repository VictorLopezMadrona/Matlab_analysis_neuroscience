function [varargout] = aw_signif_comps2(cnfg, AnyWaveData)
%% Statistical comparison between conditions
% T-test is performed on the multi-trial data detecting significant
% diferences over all the trials. Then, lfdr in applied to threshold the 
% detection.
% 
% Syntax:  
%    [I] = dyn_signif_comps(cnfg, AnyWaveData) finds statistically 
%       significant components and return their indecies.
%          .
% 
%    [I, statsum] = dyn_signif_comps(cnfg, AnyWaveData) finds statistically
%       significant components and return their indecies and number of
%       significant time samples.
% 
%    [I, statsum, t_signif] = dyn_signif_comps(cnfg, AnyWaveData) also
%       returns the statistical outcome for each component and its samples.
%
%    [I, statsum, t_signif, pval] = dyn_signif_comps(cnfg, AnyWaveData) also
%       returns the uncorrected pvalue
%
% Inputs:
%    cnfg - (Optional) structure of parameters:
%       
%       minsamples        = integer, the minimum number of consecutive 
%                           samples allowed to be detected as significant. 
%                           Default is 2, where segments of at least two 
%                           consecutive samples will be detected.
%
%                           Instead of a single value, a vector of
%                           [1 x Ncomp] can be used to apply a different
%                           value to each value.
% 
%       totalminsamples   = integer, the total minimum number of samples.
%                           Components with less significant samples will be
%                           discarded. Defaults is minsamples.
% 
%       latency           = [beg end], specify time range in seconds.
%
%       time              = (Required) when latency option is specified.
%                           A [1 x Nbr_samples] matrix of time samples.
%
%       stats             = string, if 'lfdr', find significant values
%                           based on local false discovery rate. If not 
%                           defined, the p-value would be 0.05. 
%
%       plot_lfdr         = plot lfdr stats. Def=0
%                   
%     AnyWaveData - cell of size {nb_comps x 1}. Each row is a component
%                   of [nb_trials x samples].
% 
% Outputs:
%    I - Nx1 vector with the indeces of the detected significant components
% 
%   statsum - Nx1 vector with the number of significant samples for the
%       detected components
% 
%   t_signif - [Nbr components x Nbr samples] logical matrix with the
%       statistical decision on all the input components and their samples.
%
% Other m-files required: dyn_lfdr
%
% See also: dyn_lfdr,  ttest,  bwconncomp

% Author: Christos Papageorgakis <christos.papageorgakis@gmail.com>
% Modified by: Victor Lopez Madrona <v.lopez.madrona@gmail.com>           
% License: BSD (3-clause)
% Sept. 2019; Last revision: 2-Dec-2019

% 19/08/2024: Added p-val as output
% 2/12/19: Added option to compare two conditions

nargoutchk(1,4);

% Extract data dimensions
Nev    = size(AnyWaveData,1);
nb_comps = size(AnyWaveData,2);
if nb_comps==1, nb_comps=Nev; Nev=1; end
    
nb_trials1   = size(AnyWaveData{1},1);
nb_samples  = size(AnyWaveData{1},2);

if ~exist('cnfg', 'var') || isempty(cnfg), cnfg={}; end

if ~isfield(cnfg, 'minsamples') || isempty(cnfg.minsamples)
    cnfg.minsamples = 2;
end
if ~(length(cnfg.minsamples)==1 || length(cnfg.minsamples)==nb_comp)
    error('"minsamples" length must be 1 or the number of components')
end
if ~isfield(cnfg, 'totalminsamples') || isempty(cnfg.totalminsamples)
    cnfg.totalminsamples = cnfg.minsamples;
end
if ~isfield(cnfg, 'doplot') || isempty(cnfg.doplot )
    cnfg.doplot  = false;
end
if ~isfield(cnfg, 'plot_lfdr'), cnfg.plot_lfdr  = 0; end
if ~isfield(cnfg, 'latency'), cnfg.latency = []; end
if ~isempty(cnfg.latency)
    if ~isfield(cnfg, 'time') || isempty(cnfg.time)
        error("When 'cnfg.latency' is used 'cnfg.time' parameter should be specified");
    end
end
if ~isfield(cnfg,'stats'), cnfg.stats='lfdr'; end


%% Find significant components 
% try AwSendMessage(['Computing statistics on condition: ' Cond1Name]); end

%try
    % Compute statistics for each component
    %MoyenneSurEpoques1 = nan(nb_comps, nb_samples);
    t_val = nan(nb_comps, nb_samples);
    t_signif = [];
    ttestsum = [];
    for i = 1:nb_comps
        % Mean
        %MoyenneSurEpoques1(i,:) = mean(AnyWaveData{i,1},1);
        % Type 1 error
        %pointe_erreur_type1 = std(AnyWaveData1{i,1})/sqrt(nb_trials1);
        
        % T-test
        if Nev==1
            [hhh, p_val(i,:) ,ci ,stats] = ttest(AnyWaveData{i});
        else
            [hhh, p_val(i,:) ,ci ,stats] = ttest2(AnyWaveData{1,i},AnyWaveData{2,i});
        end
        t_val(i,:)= stats.tstat;
        t_signif(i,:) = hhh;
        ttestsum(i,:) = sum( t_signif(i,:));
    end
    
    % Compute the local false discovery rate on the t-test values (not p-values)
    % over all components, to avoid the multi-test bias (repeated t-tests)
    if true
        %figure; histogram(t_val); size(t_val)
        % Compute the thresholds to apply on the t-values
        %[thr_low, thr_high] = dyn_lfdr(t_val, 0.4, 'Normal', cnfg.doplot);
        t_val_lfdr = [-t_val(:)' t_val(:)'];
        [thr_low, thr_high] = dyn_lfdr(t_val_lfdr, 0.4, 'Normal', cnfg.plot_lfdr);
        %[thr_low, thr_high] = lfdr(t_val, 0.4, 'Normal', cnfg.doplot);
        % Get the exterior points of the thresholds   <- .|.*.|. ->
        if strcmp(cnfg.stats,'lfdr')
            t_signif = ( t_val <= thr_low | t_val >= thr_high );
        else
            t_signif = p_val<0.05;
        end
        % The number of significant samples for each component
        ttestsum = sum(t_signif, 2);
        %nonzeros(ttestsum)
    end
    
    
    % Find the most statistically important components, discarding 
    % components that does not have sufficient samples per connected segment
    if cnfg.minsamples > 0
        for curr_idx = 1:nb_comps
            if length(cnfg.minsamples)>1
                minsamples=cnfg.minsamples(curr_idx);
            else
                minsamples=cnfg.minsamples;
            end
            % remove values outside the desired latency
            if ~isempty(cnfg.latency) && ~isempty(cnfg.time)
                curr_t_signif = t_signif(curr_idx,:); % nnz(curr_t_signif)
                curr_t_signif(cnfg.time{1}<cnfg.latency(1)) = 0;
                curr_t_signif(cnfg.time{1}>cnfg.latency(2)) = 0;
                t_signif(curr_idx,:) = curr_t_signif;
            end
            
            % Find the connected components in the current component
            curr_cc = bwconncomp_vector(t_signif(curr_idx,:));
            %curr_cc = bwconncomp(t_signif(curr_idx,:), 4);
            % Discard conn comps with less that minsamples
            for curr_cc_idx = 1:curr_cc.NumObjects
                if length(curr_cc.PixelIdxList{curr_cc_idx}) < minsamples
                    t_signif(curr_idx, curr_cc.PixelIdxList{curr_cc_idx}) = 0;
                end
            end
        end
        % The number of significant samples for each component
        ttestsum = sum(t_signif, 2);
    end

    
    % Get the most statistically important components, removing components
    % with low number of samples
    if true
        % Remove components without significant samples (zero counts)
        I = find(ttestsum >= cnfg.totalminsamples);     %I acts on AnyWaveData shape and similar
        ttestsum = ttestsum(I);
        
        % Determine significance in descending order of significant samples
        [~,Isorted] = sort(ttestsum,'descend');
        % Put them in the sorted order
        I = I(Isorted);
        ttestsum = ttestsum(Isorted);
        
        %     % Define the maximun number of subplots
        %     N = 20;
        %     % When only a few data exist, lower the maxumun number
        %     if N>length(ttestsum), N = length(ttestsum); end
        %     if nnz(ttestsum) < N; N = nnz(ttestsum); end
        %     % Keep only N components
        %     I = I(1:N);  % display index, samples:  [I, ttestsum(I)]
        %good_comp = I;
        % Select the SEEG components by the maximum amplitude
        if false
            % group the channels by electrode (for SEEG)
            % example: Hp6, Hp7, Hp8, Hp12 finding the 'Hp' within channels
            [~, group_indices] = findStrAlphabeticGroups(cnfg.labels(I));
            % get the component with the maximum amplitude
            I_selected = maxAmplComps(MoyenneSurEpoques1(I,:), group_indices);
            % update with the selected components only
            I=I(I_selected);
        end
        
    end
    
%catch ME
%    errordlg([ME.message ' - line: ' num2str([ME.stack(:).line]) ]);
%end

if nargout > 1
    varargout = {I, ttestsum, t_signif, p_val};
    varargout = varargout(1:nargout);
end



end