function src_out = ica_exp_var(src_in,data)

% Computes the variance that is explained by each independent component.
% It also reorders the IC matrix in function of the variance.
% AND it converts the matrix to do the loadings "positive".
% 
% USE: src_out = ica_exp_var(src_in,data);
%
% INPUT:
%   src_in: result from ft_componentanalysis
%   data:   it can be a matrix with the original data or a fieldtrip struct
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Sep. 2020; Last revision: 12-Jan-2021

if isstruct(data)
    data = data.trial{1};
end
IC = src_in.unmixing*data;
Nic = size(IC,1);

exp_var=zeros(Nic,1);
for ic=1:Nic
    data_recons = src_in.topo(:,ic)*IC(ic,:);
    exp_var(ic) = 1 - (var(data(:)) - var(data_recons(:)) ) / var(data(:));
end

[~,p]=sort(exp_var); p=p(end:-1:1);
src_out = src_in;
src_out.exp_var = exp_var(p);
src_out.topo = src_out.topo(:,p);
src_out.unmixing = src_out.unmixing(p,:);
if isfield(src_out,'label')
    src_out.label = src_out.label(p);
end

Mmax=max(src_out.topo);
Mmin=min(src_out.topo);
p=(Mmax>abs(Mmin))-(Mmax<abs(Mmin));
src_out.topo = src_out.topo.*p;
src_out.unmixing = src_out.unmixing.*p';

