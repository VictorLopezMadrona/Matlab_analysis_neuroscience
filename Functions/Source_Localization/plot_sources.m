function plot_sources(cnfg, functional, anatomical)
%% Plots functional and anatomical data while highligting locations of interest
% It highlights locations (with spheres) on the anatomical data while
% allowing to overlay the functional data as a transparent layer on top.
%
% 
% roi_loc    the location(s) of interest, to highlight on the anatomical data.
%
% Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  
%    plot_sources(cnfg, functional, anatomical) plots functionanl and
%    anatomical data while highlighting location(s) of interest placing
%    sphere(s) at the given locations on the anatomical data.
%
% Inputs:
%    cnfg - (Optional) structure of parameters:
%       
%       loi     = Nx3 vector, coordinates of the locations of interest.
%       th      = colorbar limits [min max].      
%
%       < Any option suported by ft_sourceplot and ft_sourceinterpolate >
% 
%
%    functional - functional data, see ft_sourceplot
% 
%    anatomical - anatomical data, see ft_sourceplot
%
%
% Example: 
%    cfg            = [];
%    cfg.funparameter = 'pow_nai';
%    cfg.parameter    = 'pow_nai';
%    cfg.method       = 'slice';
%    cfg.loi          = source.pos(1024,:);
%    plot_sources(cfg, source, mri_realigned);
% 
% Other m-files required: highlight_anat, ft_sourceplot,  ft_sourceinterpolate
%
% See also: ft_sourceplot,  ft_sourceinterpolate,  highlight_anat

% Author: Christos Papageorgakis <christos.papageorgakis@gmail.com>
% License: BSD (3-clause)
% Sep. 2019; Last revision: 10-Apr-2025
% Updated by Victor Lopez-Madrona

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('cnfg', 'var') || isempty(cnfg), cnfg={}; end
if ~isfield(cnfg, 'funparameter') || isempty(cnfg.funparameter)
    error("The field cfg.funparameter should be specified");
end

if isfield(cnfg, 'loi')
    % Place sphere(s) at the loi(s) on the anatomical data
    anat = highlight_anat(cnfg, anatomical);
else
    anat = anatomical;
    disp('Location(s) of interest were not specified, see cnfg.loi option.');
end

    
%%% Plot functional and highlighted anatomical data
% cfg            = [];
% cfg.parameter = 'pow';
src_interp  = ft_sourceinterpolate(cnfg, functional , anat);

% cfg              = [];
% cfg.method       = 'ortho';
% cfg.funparameter = 'pow';
ft_sourceplot(cnfg, src_interp);

% Thresholding colorbar
if isfield(cnfg, 'th') && strcmp(cnfg.method,'ortho')
    fig = gcf;                        % Get current figure
    axesList = findall(fig, 'Type', 'axes');  % Get all axes in figure
    for i = 1:length(axesList)
        ax = axesList(i);
        axes(ax);                    % Make it current (optional)
        caxis(ax, cnfg.th);            % Set your desired clim here
    end
    fig = gcf;                                      % Get current figure
    cb = findall(fig, 'Type', 'ColorBar');          % Find the colorbar handle
    if ~isempty(cb)
        set(cb, 'Limits', cnfg.th);
    end
end


end


