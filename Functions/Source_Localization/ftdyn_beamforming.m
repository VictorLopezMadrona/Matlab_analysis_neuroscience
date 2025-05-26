function [ftsource,source] = ftdyn_beamforming(cnfg,ftdata)

%% From ftdata and anatomy, obtain the time-course at each source location using beamforming
%
% Syntax:
%    [ftsource] = ft_beamforming(cnfg,ftdata)
%
% Inputs:
%    cfg - Structure of parameters:
%       fs_mri   - Path and file with MRI from freesurfer
%       spm_gm   - Gray Matter obtained from SPM 
%       grid_res - Grid resultion. Def: 10mm
%
%       dosave   - logical. True/false save/not save the figures (not the data)
%       outpath  - string. Path to save the figures if dosave=true
%       doplot   - logical. Plot the figure with the results.
%       infosave - string to include in the saved filed
%
%   ftdata - Data in FieldTrip format.
%
% Outputs:
%    ftsource - time source at each source location (FieldTrip data)
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug. 2025; Last revision: 11-Apr-2025
%

warning('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
warning('It is only working for VBMEG data, with sensors already algined with MRI')
warning('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Initialization

if ~isfield(cnfg,'grid_res'), cnfg.grid_res = 10; end

if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if ~isfield(cnfg,'outpath') && cnfg.dosave
    error('Outpath has not been specified to save the results'), end

if cnfg.dosave
    if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
    if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
end

if cnfg.dosave
    warning('Only figures will be saved. If interested in beamforming result, save it manually'), end

%%

% Get sensor position from ftdata
grad = ftdata.grad;

% Read FreeSurfer with FieldTrip
mri = ft_read_mri(cnfg.fs_mri);

%%% If VBMEG data and already aligned, define this coordsys:
mri.coordsys = 'acpc'; %Already aligned using SPM
ftdata.grad.coorsdsys = 'acpc';
grad.coorsdsys = 'acpc';

% Compute brain mesh
cfg = [];
cfg.output = 'brain';
segmented = ft_volumesegment(cfg, mri);

% Coregistration - No need in OSE dataset
%cfg = [];
%cfg.method = 'interactive';  % or 'fiducial' if you have coordinates
%cfg.headshape = 'your_meg_headshape.pos';  % FROM MEG data
%cfg.coordsys = 'freesurfer';  % because MRI is in FreeSurfer space
%mri_aligned = ft_volumerealign(cfg, mri);

% Headmodel
cfg = [];
cfg.method = 'singleshell';
headmodel = ft_prepare_headmodel(cfg, segmented);

%----------------------------------------------------------------------------------------------------------
% Plot headmodel
%----------------------------------------------------------------------------------------------------------
if cnfg.doplot || cnfg.dosave
    h=figure;
    ft_plot_sens(grad, 'unit', 'mm');
    ft_plot_headmodel(headmodel, 'facecolor', 'cortex', 'unit', 'mm');
    if cnfg.dosave
        savefig(h,[cnfg.outpath 'headmodel' cnfg.infosave ])
        saveas(h,[cnfg.outpath 'headmodel' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h), end
end
%----------------------------------------------------------------------------------------------------------

% Discretize gray matter from SPM (comptued in OSE dataset)
mri_gm = ft_read_mri(cnfg.spm_gm);
mri_gm.coordsys = 'acpc'; %Already aligned using SPM
cfg = [];
cfg.grid.resolution = cnfg.grid_res;      % 10 mm grid spacing
cfg.grid.unit = 'mm';
cfg.mri = mri_gm;
%cfg.inwardshift = -1.5;        % ensures the grid stays inside the brain surface
%cfg.tight = 'yes';
sourcemodel = ft_prepare_sourcemodel(cfg);

% Leadfield for the Nolte headmodel, created using FieldTrip
cfg                = [];
cfg.grad           = grad;
cfg.headmodel      = headmodel; 
cfg.sourcemodel    = sourcemodel; 
leadfield   = ft_prepare_leadfield(cfg);
grid = leadfield;

%----------------------------------------------------------------------------------------------------------
% compute and plot the amplitudes of the leadfields
%----------------------------------------------------------------------------------------------------------
if cnfg.doplot || cnfg.dosave
    grid_aux = {};
    grid_aux{1} = leadfield;
    
    ampl = {};
    i=1;
    ampl{i} = nan(grid_aux{i}.dim);
    for k=find(grid_aux{i}.inside(:)')
        ampl{i}(k) = sqrt(sum(grid_aux{i}.leadfield{k}(:).^2));
    end
    
    % interpolating the data to the mri for plotting
    sourceinterp = {};
    cfg             = [];
    cfg.parameter   = 'ampl';
    source          = grid_aux{i};
    source.ampl     = ampl{i};
    sourceinterp{i} = ft_sourceinterpolate(cfg, source, mri);
    
    % plotting the amplitudes
    h=figure;
    cfg               = [];
    cfg.funparameter  = 'ampl';
    cfg.method        = 'slice';
    ft_sourceplot(cfg, sourceinterp{1});
    title('Amplitude of leadfield')
    
    if cnfg.dosave
        %savefig(h,[cnfg.outpath 'leadfield' cnfg.infosave ])
        saveas(h,[cnfg.outpath 'leadfield' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h), end
end
%----------------------------------------------------------------------------------------------------------

%%% Beamforming

%%% Compute covariance matrix - raw data
cfg                  = [];
cfg.covariance       = 'yes';   % required for the lcmv beamformer: Cov matrix
cfg.covariancewindow = 'all';   % use all signal
cfg.vartrllength     = 0;
rawtimelock          = ft_timelockanalysis(cfg, ftdata);

%%% Compute beamforming
cfg        = [];
cfg.grid = grid; % Use precomputed leadfield
cfg.headmodel = headmodel;
cfg.method = 'lcmv';
cfg.lcmv.reducerank = 'no'; %default MEG = 2
cfg.lcmv.lambda       = '10%'; % Velu parametre 10%
cfg.lcmv.fixedori     = 'yes'; %only one direction instead of three
cfg.lcmv.keepfilter   = 'yes';
cfg.lcmv.keepmom      = 'yes';
cfg.lcmv.projectnoise = 'yes';
source = ft_sourceanalysis(cfg, rawtimelock);

%%% Neural Activity Index (NAI)
source.avg.pow_nai = source.avg.pow./source.avg.noise;

%----------------------------------------------------------------------------------------------------------
% Plot beamforming Neural Activity Index (NAI)
%----------------------------------------------------------------------------------------------------------
if cnfg.doplot || cnfg.dosave
    cfg = [];
    cfg.funparameter = 'pow_nai';
    cfg.parameter = 'pow_nai';
    cfg.method = 'ortho'; % surface - ortho - slice
    cfg.th = [max(source.avg.pow_nai)*0.5 max(source.avg.pow_nai)];
    h=figure; 
    plot_sources(cfg,source,mri)
    if cnfg.dosave
        %savefig(h,[cnfg.outpath 'beamforming_nai' cnfg.infosave ])
        saveas(h,[cnfg.outpath 'beamforming_nai' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h), end
end
%----------------------------------------------------------------------------------------------------------

% Time-course at each location
ftsource = conv_source_to_raw(source);





