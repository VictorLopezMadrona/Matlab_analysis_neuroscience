function ftdata_noERP = removeERP_cont(cnfg,ftdata)

%% Function to substract the ERP from each trial, but in a continuous data
% 1- Obtain trials from continuous data
% 2- Use the ERP as a regressor for each trial.
% 3- Subtract the regressed component to get the residual signal.
% 4- Reconstruct continuous data
%
% ftdata_noERP = removeERP_cont(cnfg,ftdata)
%
% INPUT: ftdata with trials
%   cfg.stimdef    - [start end] Trial in seconds. Ex [-0.2 0.5]
%       time_trial - vector with each trial onset in seconds
%
% OUTPUT:
%   ftdata_noERP: the same as ftdata without ERPs
%   amp_ERP: amplitude of the ERPs for each trial (obtained from the regression)
%
% Victor J. Lopez-Madrona
% 11/03/2025

%%

% 1- Obtain trials from continuous data
stimdef    = cnfg.stimdef;
time_trial = cnfg.time_trial;
Ntrial  = length(time_trial);
ft_time = ftdata.time{1};
Fs      = ftdata.fsample;
Nch = size(ftdata.trial{1},1);

onset_samples = zeros(1,Ntrial);
for tri = 1:Ntrial
    [~, onset_samples(tri)] = min(abs(ft_time - time_trial(tri)));  % index of minimum difference
end
pre_samples   = round(stimdef(1)*Fs); % Samples before onset in samples
post_samples  = round(stimdef(2)*Fs); % Samples after onset in samples
% Define trials
trl_ini = onset_samples + pre_samples;  % Start sample
trl_end = onset_samples + post_samples; % End sample
Nsamples = trl_end(1)-trl_ini(1)+1;
data_trial = zeros(Nch,Nsamples,Ntrial);
for tri=1:Ntrial
    data_trial(:,:,tri) = ftdata.trial{1}(:,trl_ini(tri):trl_end(tri));
end

% 2- Compute the averaged ERP
ERP_avg = mean(data_trial,3);

% 3- Substract the ERP from each trial and channel
data_trial_noERP = zeros(Nch,Nsamples,Ntrial);
for tri=1:Ntrial
    for chi=1:Nch
        % Least squares
        A = ERP_avg(chi,:); A=A(:);
        B = data_trial(chi,:,tri); B=B(:);
        % Perform linear regression to estimate the contribution of A to B
        amp_ERP = (A' * A) \ (A' * B);  % Least squares solution 
        % Calculate the predicted ERP contribution
        erp_contribution = A * amp_ERP;
        % Subtract the ERP contribution to get the residual signal
        data_trial_noERP(chi,:,tri) = B - erp_contribution;
    end
end

% 4- Reconstruct continuous data
ftdata_noERP = ftdata;
for tri=1:Ntrial
    ftdata_noERP.trial{1}(:,trl_ini(tri):trl_end(tri)) = data_trial_noERP(:,:,tri);
end






