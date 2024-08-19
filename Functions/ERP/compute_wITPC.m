function [wITPCz,wITPC,pval] = compute_wITPC(k,A,Ns)

%% Computed the weighted ITPC with stats using permutation test
%
% USE:
%   [wITPCz,wITPC,pval] = compute_wITPC(phase,weight,Nperm);
%
% INPUT:
%   phase = vector with the phases
%   weight = vector with the weights associated to each phase
%   Nperm = Number of permutations for statistical test (Nperm = 1000 is ok)
%
% OUTPUT:
%   wITPCz = z-scored weighted ITPC result. Corrected using the permutations
%   wITPC = weighted ITPC result. Not very informative
%   pval
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Aug. 2024; Last revision: 16-Aug-2024

wITPC = abs(mean(A.*exp(1i*k)));

wITPCs = zeros(1,Ns+1); % Number of permutations
for si=1:Ns
    As = A(randperm(length(A)));
    wITPCs(si) = abs(mean(As.*exp(1i*k)));
end
wITPCs(Ns) = wITPC; %I add the observed value into the dist
pval = sum(wITPCs>=wITPC)/Ns;
wITPCz = (wITPC-mean(wITPCs)) / std(wITPCs);




