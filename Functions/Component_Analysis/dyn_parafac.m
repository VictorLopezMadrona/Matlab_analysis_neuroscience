function [F,T,W] = dyn_parafac(cnfg,data3d)

%% From 3D data compute parafac. Requires N-Way toolbox
% Ideally for frequency x time x space (sensors) data.
% Only plots CORCONDIA results. Other plottings must be done after.
%
% Syntax:
%    [T,F,W] = dyn_parafac(cnfg,data3d)
%
% Inputs:
%    cfg - Structure of parameters:
%
%       K          - [Optional if Kmax] Number of components
%       Kmax       - [Optional if K] Find best K using corcondia
%       stats      - method for stats: Def 'std'
%       alpha      - pval for std stats. Def= 0.01
%       corcond_th - Threshold for corcondia in % (Def 80)
%       cons       - Constraints for each dimension. Def [0 0 2]: 
%                       0- no constraints
%                       2- nonnegative constraint
%                   See 'parafac.m' for more cases
%
%       dosave   - logical. True/false save/not save the results
%       outpath  - string. Path to save the results if dosave=true
%       doplot   - logical. Plot the figure with the results.
%       infosave - string to include in the saved filed
%
%   data3d - time-frequency-space data
%
% Outputs:
%    T   - Temporal distribution of each mode (component)
%    F   - Frequency distribution of each mode (component)
%    W   - Spatial distribution of each mode (component)
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug. 2025; Last revision: 15-Apr-2025

% Changes log:

%% Initialization

if ~isfield(cnfg,'K') && ~isfield(cnfg,'Kmax')
    error('It is mandatory to define K or Kmax')
end
if isfield(cnfg,'K') && isfield(cnfg,'Kmax')
    warning('Both K and Kmax were defined')
    warning('Using K instead of obtaining optimal K')
    cnfg = rmfield(cnfg, 'Kmax');
end
if ~isfield(cnfg,'cons'), cnfg.cons = [0 0 2]; end
if ~isfield(cnfg,'corcond_th'), cnfg.corcond_th = 80; end
if ~isfield(cnfg,'stats'), cnfg.stats = 'std'; end
if ~isfield(cnfg,'alpha'), cnfg.alpha = 0.01; end
alpha = norminv(1-cnfg.alpha);

if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if ~isfield(cnfg,'outpath') && cnfg.dosave
    error('Outpath has not been specified to save the results'), end

if cnfg.dosave
    if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
    if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
end

%% PARAFAC

%%% Find best K based on CORCONDIA
if isfield(cnfg,'Kmax')
    corcond = zeros(1,cnfg.Kmax);
    for ki=1:cnfg.Kmax
        [~,~,~,corcond(ki)] = parafac(data3d,ki,[],cnfg.cons);
    end
    
    % Define optimal K
    cnfg.K = find(corcond >= cnfg.corcond_th,1,'last');
    
    % Plot corcondia
    h=figure;
    hold on,
    plot(corcond,'LineWidth',1.5)
    plot([0 cnfg.Kmax+1],[0 0],'k','LineWidth',1)
    plot([0 cnfg.Kmax+1],[cnfg.corcond_th cnfg.corcond_th],'--k','LineWidth',1)
    for ki=1:cnfg.Kmax
        plot([ki ki],[0 100],'k','LineWidth',0.5)
    end
    plot(cnfg.K,corcond(cnfg.K),'rx','MarkerSize',15,'LineWidth',2)
    axis([0.5 cnfg.Kmax+0.5 0 100])
    xticks(1:cnfg.Kmax)
    xlabel('Number of modes (k)')
    ylabel('CORCONDIA value')
    
    if cnfg.dosave
        %savefig(h,[cnfg.outpath 'Corcondia' cnfg.infosave ])
        saveas(h,[cnfg.outpath 'Corcondia' cnfg.infosave '.png'])
    end
    if ~cnfg.doplot, close(h), end
end

%%% Perform PARAFAC

Factors = parafac(data3d,cnfg.K,[],cnfg.cons);
[F,T,W] = fac2let(Factors);


