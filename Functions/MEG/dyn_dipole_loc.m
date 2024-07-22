function [loc,channel,MaxGoF,GoF_SEEG]=dyn_dipole_loc(cnfg)

%% Dipole localization of ICA topographies
% 
% Syntax:  
%    [loc,channel]=dyn_dipole_loc(cfg);
%
% Inputs:
%    cfg - Structure of parameters:
%       
%       mri_realigned = path with mri_realigned
%       leadfield = path with leadfield
%       ica = path the ICA topography
%       mode = 'single', '2sym'
%       channel = ICs to analyze. Def=all.
%       SEEGpos = location to compute GOF there
%       doplot = True/False (Def=False);
%
% Outputs:
%    loc     = Location of each dipole
%    channel = IC associated to each location        
%
% See also: write_graphcorr_oldnew

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% May. 2020; Last revision: 6-Apr-2021

if ~isfield(cnfg,'doplot'), cnfg.doplot=false; end
    
grid=[];
load(cnfg.leadfield)
load(cnfg.ica)
load(cnfg.mri_realigned)

if ~isfield(cnfg,'channel')
    channel = 1:size(mixing,2);
else
    channel = cnfg.channel;
end

comp.topo = mixing;

[new_leadfield, new_labels] = dyn_reordering_leadfield_to_current_channels(labels, grid);
grid.leadfield = new_leadfield;
grid.label = new_labels;

label_common=ismember(labels,new_labels);
comp.topo=comp.topo(label_common,:);

%%%%% ===== localisation avec intervalle de confiance ==== %%%%%
for i=1:length(channel)
    compi=channel(i);
    CompOI = compi; % composante d'interet
    cur_data = comp.topo;
    L = grid.leadfield;
    ssq_total=sum(sum(cur_data(:, CompOI).^2));
    idx_inside = find(grid.inside);
    n_inside = length(idx_inside);
    ssq1=zeros(n_inside,1);
    
    topo=cur_data(:,CompOI);
    
    if strcmp(cnfg.mode,'2sym')
        for j=1:length(idx_inside)
            curj=idx_inside(j);
            curj_op=find(grid.pos(:,1)==grid.pos(curj,1) & ...
                grid.pos(:,3)==grid.pos(curj,3) & grid.pos(:,2)==-grid.pos(curj,2));
            curL1=[cell2mat(L(curj)) cell2mat(L(curj_op))];
            %curL1=[cell2mat(L(curj))];
            S1=pinv(curL1)*topo;
            model1=curL1*S1;
            ssq1(j)=sum((topo(:)-model1(:)).^2);
        end
    else
        for j=1:length(idx_inside)
            curj=idx_inside(j);
            curL1=[cell2mat(L(curj))];
            S1=pinv(curL1)*topo;
            model1=curL1*S1;
            ssq1(j)=sum((topo(:)-model1(:)).^2);
        end
    end
    
    % meilleure poistion du single dipole fit
    [~,I]=min(ssq1(:));
    PosScan1dipole = idx_inside((I));
    
    % gof :
    gof1dipole=1-ssq1/ssq_total;
    MaxGoF(i) = max(gof1dipole);
    
    loc(i,:)=grid.pos(PosScan1dipole,:);
    
    if isfield(cnfg,'SEEGpos')
    for si=1:size(cnfg.SEEGpos,1)
        dist=sqrt(sum( (grid.pos(idx_inside)-cnfg.SEEGpos(si,:))'.^2 ));
        [~,p]=min(dist);
        GoF_SEEG(si,i)=gof1dipole(p);
    end
    else
        GoF_SEEG=0;
    end
       
    % ======= display 2 : Intervalle de confiance
    if cnfg.doplot
        ICrange = [max(0, MaxGoF(i) - (1-MaxGoF(i))) , MaxGoF(i)];
        if ICrange(1)<0.8 && ICrange(2)>0.8
            ICrange = [0.8 , MaxGoF(i)];
        end
        
        ToInterpolate = [];
        ToInterpolate.inside = grid.inside;
        ToInterpolate.pos = grid.pos;
        ToInterpolate.dim = grid.dim;
        pow = zeros(size(grid.inside,1),1);
        pow(idx_inside,:) = gof1dipole.*(gof1dipole>ICrange(1,1));
        ToInterpolate.avg.pow = pow;
        cfg            = [];
        cfg.parameter = 'pow';
        cfg.funcolormap = 'parula';
        Interpolated1 = ft_sourceinterpolate(cfg,ToInterpolate, mri_realigned);
        
        cfg = [];
        cfg.method = 'ortho';
        cfg.location = grid.pos(PosScan1dipole,:);
        cfg.interactive = 'yes';
        cfg.crosshair = 'yes';
        cfg.funparameter = 'pow';
        cfg.maskparameter = cfg.funparameter;
        if min(ICrange)>0
            cfg.opacitylim = min(ICrange)*[1 1.0000001];
        end
        cfg.axis = 'on';
        cfg.funcolormap = 'parula';
        cfg.funcolorlim = ICrange;
        %cfg.visible = 'off';
        ft_sourceplot(cfg, Interpolated1);
    end
    
end

    

