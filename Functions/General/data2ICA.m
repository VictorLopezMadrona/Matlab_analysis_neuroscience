function ftdata=data2ICA(cnfg,ftdata)

%% Creates a FT struct with ICA on ftdat
% 
% Syntax:  
%    ftdataIC=data2ICA(cfg,ftdata);
%
% Inputs:
%   cfg.ICname  = Path with the ICA file 
%   cfg.src_ica = struct with ica
%   ftdata = struct as the one given by ft_preprocessing.
%
% Outputs:
%   ftdataIC = FT struct changing trials and labels.
%
% See also: loadMEGdata

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Nov. 2019; Last revision: 30-Nov-2020


if ~isfield(cnfg,'ICname') && ~isfield(cnfg,'src_ica')
    error('A path with the ICA file or src_ica must be specified.')
end
if isfield(cnfg,'ICname') && ~exist(cnfg.ICname,'file')
    error("Selected ICA file doesn't exist")
end

if isfield(cnfg,'src_ica')
    src_ica = cnfg.src_ica;
else
    load(cnfg.ICname)
end

%Anywave option
if exist('unmixing','var') && exist('labels','var')
%FT option: src_ica
elseif exist('src_ica','var')
    unmixing = src_ica.unmixing;
    labels = src_ica.topolabel;
else
    error('Do not recognize ICA input information')
end

% Change the data
[~,ord]=ismember(ftdata.label,labels);

%unmixing=unmixing(:,ord);
unmixing=unmixing(:,ord(ord>0));
[nic,nch_ic]=size(unmixing);
[nch,~]=size(ftdata.trial{1});
if nch_ic~=nch
    %error('Number of channels and IC matrix must agree')
    ftdata_ch = find(ord>0);
else
    ftdata_ch = 1:nch;
end
    
%Create labels
if exist('labelsICA','var')
    ftdata.label = labelsICA';
else
    ICAlabel=cell(nic,1);
    for nic_i=1:nic
        ICAlabel{nic_i,1} = ['ICA ' num2str(nic_i)];
    end
    ftdata.label = ICAlabel;
end

%Change Channels by IC
for trl_i = 1:length(ftdata.trial)
    ftdata.trial{trl_i} = unmixing*ftdata.trial{trl_i}(ftdata_ch,:);
end


