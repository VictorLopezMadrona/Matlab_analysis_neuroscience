function source = source_pow_window(cnfg, source, src_as_raw)
%% Computes the power of a selected time window
% For the multi-trial case, the power is computed on the average across
% trials for the dara within the selected time window.
%
% Syntax:  
%    src = source_pow_window(cnfg, source) computes the power or a selected
%    time window.
% 
%    src = source_pow_window(cnfg, source, src_as_raw) computes the power 
%    of a selected time window using precomputed src_as_raw
% 
% Inputs:
%    cnfg - (Optional) structure of parameters:
% 
%       latency = [begin end] in seconds, the time window of interest
%       
%       trl     = Nx3 matrix with the trial definition, see FT_DEFINETRIAL,
%                 required when the input source not in in a raw like
%                 structure.
% 
%    source - functional data, see ft_sourceplot
% 
%    src_as_raw - (Optional) functional data converted to raw, see 
%                 conv_source_to_raw
%
% Outputs:
%    src - A source structure with added fields: 
% 
%          pow_win        the power of the data within the time window
% 
%          pow_win_nai    the neural activity index (power) of the data 
%                         within the time window
%
% Example: 
%    cfg = [];
%    cfg.trl = trl;
%    cfg.latency = [0.1, 0.12];
%    source_pow_window(cfg, source)
%
% Other m-files required: conv_source_to_raw, ft_redefinetrial
%
% See also: conv_source_to_raw,  ft_redefinetrial,  ft_timelockanalysis

% Author: Christos Papageorgakis <christos.papageorgakis@gmail.com>
% License: BSD (3-clause)
% Oct. 2019; Last revision: 16-Oct-2019

% The power is computed as the power of the average over trials.

if ~exist('cnfg', 'var') || isempty(cnfg), cnfg={}; end
if ~isfield(cnfg, 'latency'), cnfg.latency = []; end

if ~ft_datatype(source, 'source')
    error('Input "source" should be of a source structure.')
end

%%% Source -> Raw -> Epochs
if ~exist('src_as_raw', 'var') || ft_datatype(source, 'raw')
    % Source -> Raw
    disp('Converting source to raw structure...');
    src_as_raw = conv_source_to_raw(source);
else
    disp('Using precomputed "src_as_raw"...');
end

if isfield(cnfg, 'trl') && length(src_as_raw.trial)<2
    disp('Converting raw structure to epochs...');
    % Raw -> Epochs
    cfg = [];
    cfg.trl = cnfg.trl;
    src_as_epochs = ft_redefinetrial(cfg, src_as_raw);
else
    src_as_epochs = src_as_raw;
end

clear src_as_raw;


%%% Select the data within the time window
cfg      = [];
cfg.latency = cnfg.latency;
src_as_epochs_win_timelock = ft_timelockanalysis(cfg, src_as_epochs);

clear src_as_epochs;

%%% Compute the power of the selected time window
pow_win = sum(src_as_epochs_win_timelock.avg .^2, 2);

%%% Store the power of the selected time window
source.avg.pow_win = source.avg.pow;
source.avg.pow_win(source.inside) = pow_win;

%%% Compute the NAI power of the selected time window
source.avg.pow_win_nai = source.avg.pow_win;

%%% Store the NAI power of the selected time window
source.avg.pow_win_nai = source.avg.pow_win_nai ./ source.avg.noise;


end