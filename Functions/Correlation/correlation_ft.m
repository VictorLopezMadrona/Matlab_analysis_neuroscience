function [maxcor,del,cor,lag,THpos,THneg]=correlation_ft(cnfg,ftdata)

%% Compute cross-correlation between selected channels in ft structure
%
%  Remember: a positive delay means that the second signal preceeds the first 
%
% USE:
%   [cor,del]=correlation_ft(cfg,ftdata)
%
% PARAMETERS (cfg.): 
%   
%   latency     - [t_start t_end]
%   window      - window length in time (s) [Def=5]
%   step        - step between windows in time (s) [Def=window]
%   maxlag      - max lag for xcorr (improve computational time) [Def=window-1]
%   resample    - Fs to resample data after filtering
%   Nsurro      - Number of surrogates to statistical significance (Def=100)
%   doplot      - true/false (Def: true)
%   dosave      - true/false (Def: false)
%   outpath     - string with path 
%   infosave    - string with filename
%
% OUTPUT: Modify this
%   maxcor - max value of xcor (pos or neg) for each pair of signals
%   del    - delay associated to maxcor
%   cor    - dynamics of the xcor within maxlag averaged across windows
%   lag    - lag associated to each position in 'cor'
%   THpos  - threshold of significance for positive values
%   THpos  - threshold of significance for negative values
%
% See also: 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Feb 24; Last revision: 18-Feb-2025

%% TO DO:
% Include the option to plot the cross-correlation in time for each pair
% Include the option to plot the distribution of delays for each time
% across windows
%

%% Initial parameters:

if ~isfield(cnfg,'latency')
    cnfg.latency=[ftdata.time{1}(1) ftdata.time{1}(end)]; end
if ~isfield(cnfg,'window'), cnfg.window = 5; end
if ~isfield(cnfg,'step'), cnfg.step = cnfg.window; end
if ~isfield(cnfg,'maxlag')
    cnfg.maxlag = (cnfg.window*ftdata.fsample)-1;
    warning('It is recommended to define maxlag to reduce compational cost')
end
if ~isfield(cnfg,'Nsurro'), cnfg.Nsurro = 100; end
if ~isfield(cnfg,'pval'), cnfg.pval =   0.05; end
if ~isfield(cnfg,'Nsurro'), cnfg.Nsurro = 0; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if isfield(cnfg,'outpath')
    if ~strcmp(cnfg.outpath(end),'\')
        cnfg.outpath = [cnfg.outpath '\'];
    end
end
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end
cnfg.Fs = ftdata.fsample;

%% Resample

Fs = ftdata.fsample;
if isfield(cnfg,'resample'), Fs = cnfg.resample; end
data = ftdata.trial{1};
if isfield(cnfg,'resample')
    for chi=1:size(data,1)
        data_aux(chi,:) = resample(data(chi,:),cnfg.resample,round(cnfg.Fs));
    end
    data=data_aux;
    clear data_aux
    Fs = ftdata.fsample;
end
%data=zscore(data')';

%% COMPUTE XCORR

w    = round(cnfg.window*Fs);
step = round(cnfg.step*Fs);
Nw   = floor((size(data,2)-w)/step)+1;
Nch  = size(data,1);

% Compute the xcorr at each window. Then average the windows and find the
% max and its associated delay.
% Compute surrogates to test if the corr value is significant or not. If
% so, we can trust the delay. If not, no.
cor = zeros(Nch,Nch,cnfg.maxlag*2+1); % this is my result
del = zeros(Nch,Nch); % this is my result
maxcor = zeros(Nch,Nch); % this is my result
for chi=1:Nch-1
    for chj=chi+1:Nch
        C = zeros(Nw,cnfg.maxlag*2+1);
        for wi=1:Nw
            x = data(chi,(wi-1)*step+1:(wi-1)*step+w);
            y = data(chj,(wi-1)*step+1:(wi-1)*step+w);
            [C(wi,:),lag] = xcorr(x,y,cnfg.maxlag,'coeff');
        end
        cor(chi,chj,:) = mean(C);
        [~,p] = max(abs(mean(C)));
        maxcor(chi,chj) = C(p);
        maxcor(chj,chi) = C(p);
        del(chi,chj) = lag(p)/Fs*1000; %delay in ms    
        del(chj,chi) = lag(p)/Fs*1000; %delay in ms   

%         % Fast stats. Check if the selected point is sig diff from zero
%         % across windows -- Useless
%         [~,pval_ttest] = ttest(C(:,p));
%         if pval_ttest>0.05
%             maxcor(chi,chj) = 0;
%             maxcor(chj,chi) = 0;
%             del(chi,chj) = 0;
%             del(chj,chi) = 0;
%         end
    end
end

%% SURROGATES
Nsurro = cnfg.Nsurro;
THpos = zeros(Nch,Nch); % this is my result
THneg = zeros(Nch,Nch); % this is my result
for chi=1:Nch-1
    for chj=chi+1:Nch
        Csurro = zeros(Nsurro,cnfg.maxlag*2+1);
        for si=1:Nsurro
            Caux = zeros(Nw,cnfg.maxlag*2+1);
            for wi=1:Nw
                ws=randperm(Nw); % I select two random windows
                x = data(chi,(ws(1)-1)*step+1:(ws(1)-1)*step+w);
                y = data(chj,(ws(2)-1)*step+1:(ws(2)-1)*step+w);
                [Caux(wi,:),lag] = xcorr(x,y,cnfg.maxlag,'coeff');
            end
            Csurro(si,:) = mean(Caux);
        end
        %For each surrogate, I will do the average across windows and keep
        %the highest and lowest values.
        pos_val = max(Csurro,[],2);
        neg_val = min(Csurro,[],2);
        THpos(chi,chj) = norminv(1-cnfg.pval/2,mean(pos_val),std(pos_val));
        THneg(chi,chj) = norminv(cnfg.pval/2,mean(neg_val),std(neg_val));
        THpos(chj,chi) = THpos(chi,chj);
        THneg(chj,chi) = THneg(chi,chj);
    end
end

%% Plot results

mask = (maxcor>=THpos) | (maxcor<=THneg);
maxcor_mask = maxcor.*mask;
del_mask = del.*mask;
hcor=figure;
imagesc(maxcor_mask);
axis('xy')
xticks(1:Nch);
xticklabels(ftdata.label)
yticks(1:Nch);
yticklabels(ftdata.label)
colorbar
title('Cross-correlation')

hdel=figure;
imagesc(del_mask);
axis('xy')
xticks(1:Nch);
xticklabels(ftdata.label)
yticks(1:Nch);
yticklabels(ftdata.label)
colorbar
title('Delay Cross-correlation')

if cnfg.dosave
    if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
    if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
    savefig(hcor,[cnfg.outpath 'Xcorr' cnfg.infosave])
    saveas(hcor,[cnfg.outpath 'Xcorr' cnfg.infosave '.png'])
    savefig(hdel,[cnfg.outpath 'Delay_Xcorr' cnfg.infosave])
    saveas(hdel,[cnfg.outpath 'Delay_Xcorr' cnfg.infosave '.png'])    
end
if ~cnfg.doplot
    close(hcor)
    close(hdel)
end
if cnfg.dosave
    save([cnfg.outpath 'Xcorr' cnfg.infosave],'cor','del','maxcor','lag','THpos','THneg')
end


