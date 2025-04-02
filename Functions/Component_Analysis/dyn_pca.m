
function [data_pca,eigenvectors,eigenvalues]=dyn_pca(data,ncomps)

%% Computes principal component analysis following runica.m code.
% Pipeline for ICA:
% 1- Reduce number of dimensions in 'data' with dyn_pca.
% 2- Compute ICA on 'data_pca'.
% 3- The ICA matrices for all channels are:
%       > unmixing_all_channels = unmixing_pca * eigenvectors(:,1:ncomps)';
%       > mixing_all_channels   = pinv(unmixing_all_channels)
% 
% Syntax:  
%   [data_pca,eigenvectors,eigenvalues]=dyn_pca(data,ncomps);
%
% Inputs:
%   data
%   ncomps
%
% Outputs:
%   data_pca
%   eigenvectors
%   eigenvalues
%
% See also: runica.m

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% May 2020; Last revision: 5-May-2020


PCdat2 = data';                    % transpose data
[PCn,~]=size(PCdat2);                  % now p chans,n time points
PCdat2=PCdat2/PCn;
PCout=data*PCdat2;
clear PCdat2;

[PCV,PCD] = eig(PCout);                  % get eigenvectors/eigenvalues
[PCeigenval,PCindex] = sort(diag(PCD));
PCindex=rot90(rot90(PCindex));
PCEigenValues=rot90(rot90(PCeigenval))';
PCEigenVectors=PCV(:,PCindex);
%PCCompressed = PCEigenVectors(:,1:ncomps)'*data;
data_pca = PCEigenVectors(:,1:ncomps)'*data;

eigenvectors=PCEigenVectors;
eigenvalues=PCEigenValues;

clear PCn PCp PCout PCV PCD PCeigenval PCindex PCEigenValues PCEigenVectors data



