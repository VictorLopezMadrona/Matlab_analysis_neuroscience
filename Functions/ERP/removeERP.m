function [ftdata_noERP, amp_ERP] = removeERP(ftdata)

%% Function to substract the ERP from each trial using regression.
% 1- Use the ERP as a regressor for each trial.
% 2- Subtract the regressed component to get the residual signal.
%
% [ftdata_noERP, amp_ERP] = removeERP(ftdata)
%
% INPUT: ftdata with trials
%
% OUTPUT:
%   ftdata_noERP: the same as ftdata without ERPs
%   amp_ERP: amplitude of the ERPs for each trial (obtained from the regression)
%
% Victor J. Lopez-Madrona
% 11/03/2025

%% To do
% Include the option to work with data matrix instead of ft structure

% 1- Compute the averaged ERP
cfg=[];
tl_data = ft_timelockanalysis(cfg, ftdata);

ntr = size(ftdata.trial,2);
nch = size(ftdata.label,1);
nsamples = size(ftdata.trial{1},2);

data_noERP = zeros(nch,ntr*nsamples);
ftdata_noERP = ftdata;
amp_ERP = zeros(nch,ntr);

% 2- Substract the ERP from each trial and channel
for tri=1:ntr
    for chi=1:nch
        lm = fitlm(tl_data.avg(chi,:),ftdata.trial{tri}(chi,:));
        %In matrix
        data_noERP(chi,(tri-1)*nsamples+1:tri*nsamples) = lm.Residuals.Raw-mean(lm.Residuals.Raw);
        %In ft structure
        ftdata_noERP.trial{tri}(chi,:) = lm.Residuals.Raw-mean(lm.Residuals.Raw);
        %Relative amplitude of ERP
        amp_ERP(chi,tri) = lm.Coefficients.Estimate(2,1);
    end
end
