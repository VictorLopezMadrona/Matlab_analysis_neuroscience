function [thr_low, thr_high] = dyn_lfdr(xi, alpha, distname, doplot)
%% Local false discovery rates
%
% Efron, B. (2005). Local false discovery rates.
%
% Syntax:  
%    [thr_low, thr_high] = lfdr(xi, alpha, type, doplot) computes the local
%    false discovery rates of the data x and returns the minimum and
%    maximum threshold of the input data (threshold values correspond
%    to the same units as the input data).
%
% Inputs:
%    xi - Sample data, specified as a vector, matrix, or multidimensional 
%       array.
% 
%    alpha - Significance level of the hypothesis test (default is 0.05).
%
%    distname - Distribution name, specified as character vector or string
%       scalar (default is 'Normal').
% 
%    doplot - Flag to plot the distribution fit and thresholds, specified
%       as logical (defalt is 0).
%
% Outputs:
%    thr_low - The lower threshold value (return -Inf when not found).
%
%    thr_high - The higher threshold value (return Inf when not found).
% 
% Overview:
%    - Removes the outliers in the given data (=> cleaned data)
%    - Fits on 'cleaned data' the f0 distribution
%    - Computes the thresholds for the given 'alpha' precision
%
% Example: 
%     mu0 = 5; % mean of H0
%     mu1 = 8; % mean of H1
%     sigma0=0.5; % scale of H0
%     sigma1=1; % scale of H1
%     N0 = 5000; % number of sample under H0
%     N1 = 300; % number of sample under H1
%     H0 = randn(N0,1)*sigma0 + mu0; % draw N0 samples from normally distributed H0 of mean and std of (mu0,sigma0)
%     H1 = randn(N1,1)*sigma1 + mu1; % draw N1 samples from normally distributed H1 of mean and std of (mu1,sigma1)
%     x = [H0;H1];
%     [thr_low, thr_high] = lfdr(x, 0.2, 'Normal',1);
%
% See also: ttest, fitdist

% Author: Marmaduke Woodman (INS)
%		  Nicolas Roehri (INS)
% License: All rights reserved. 
% Apr. 2019; Last revision: 12-Apr-2019

if nargin < 2, alpha = 0.05; end
if nargin < 3, distname = 'Normal'; end
% if nargin < 3, distname = 'GM'; end
if nargin < 4, doplot = 0; end


% find f0
xi = xi(:);
K = 1.5;
qs = [quantile(xi,1/4) - K*iqr(xi); K*iqr(xi) + quantile(xi,3/4)];% Inf*[-1 1];
% qs = [quantile(xi,0.25); quantile(xi,0.75)];
N_smooth = 3;

if doplot
    figure;
    h = histogram(xi);
else
    % Equivalent computation without the plot    
    h=struct();
    [h.Values,h.BinEdges] = histcounts(xi);
end

% h = histogram(xi, -10:0.1:10); % 1st test on MEG
x = mean([h.BinEdges(1:end-1);h.BinEdges(2:end)]);

n = h.Values;
[~, id] = max(n); % get the mode of the distribution
M = h.BinEdges(id); % get the mode of the distribution
% n(n == 1) = 0; % get rid of outliers that would lower the threshold
f = (n + 0) / sum(n);
f = filtfilt(1/N_smooth*ones(N_smooth,1),1,f);

switch distname
    
    case 'GM'
        options = statset('MaxIter',200);
        pd = fitgmdist(xi(xi < qs(end) & xi > qs(1)),2,'Options',options);

    otherwise
        pd = fitdist(xi(xi < qs(end) & xi > qs(1)),distname);
end

f0 = pdf(pd,x');
f0 = (f0' / max(f0)) * max(f);
lf = log(f0./f);

% Find the threshold values: using the mean of each histogram bin
thr_low = x(find(lf<log(alpha) & x < M, 1,'last'));
thr_high = x(find(lf<log(alpha) & x > M, 1,'first'));


% Find the threshold values: using the mean of each histogram bin, but
% using the beggining of the bin (for each side)
if isempty(thr_low)
    thr_low = -Inf;
else
    thr_low = h.BinEdges( find(h.BinEdges > thr_low, 1, 'first') );
end
if isempty(thr_high)
    thr_high = Inf;
else
    thr_high = h.BinEdges( find(h.BinEdges < thr_high, 1, 'last') );
end



% disp([thr_low;thr_high])

% RMSE = sqrt(mean((f((x < qs(end) & x > qs(1)))-f0(x < qs(end) & x > qs(1))).^2));
% NRMSE = RMSE/max(f);
% GoF = 1 - sum((f((x < qs(end) & x > qs(1)))-f0(x < qs(end) & x > qs(1))).^2)/std(f((x < qs(end) & x > qs(1)))-f0(x < qs(end) & x > qs(1)));


%% plot results
if doplot
    figure
    ax(1) = subplot(2,1,1);
    plot(x, f, 'k');
    hold on
    plot(x, f0, 'b');
    grid on
    plot([qs(1) qs(end)],[0 0] , 'rs', 'Markersize',20, 'Linewidth',2)
    plot([thr_low thr_high],[0 0] , 'bx', 'Markersize',20, 'Linewidth',2)
    title('lfdr debugger')
    hold off
%     xlim([-1 1]);
    legend('H_G','H_0','Outlier thr', 'lfdr thr')
    ax(2) = subplot(2,1,2);
    plot(x, log(f), 'k');
    hold on
    plot(x, log(f0), 'b');
    plot(x, lf, 'r');
%     ylim([-10, 0.5]);
%     xlim([-1 1]);
    grid on
    title('lfdr debugger')
    hold off
    linkaxes(ax, 'x')
    legend('H_G','H_0','lFDR')
end

end