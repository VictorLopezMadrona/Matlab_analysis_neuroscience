
function [Favg,favg]=computeGC_MVGC(X,cnfg)

% Function to compute temporal and spectral Granger Causality
%
% USE:
%   computeGC_MVGC(X,cfg)
%
% INPUT:
%   X [nch x samples]: Time-series. Each channel should be in a different row.
%   cfg - struct with parameters:
%        .label: [Optional]
%        .window: length of the sliding window in samples;
%        .step: step of the sliding window in samples. 
%        .p: Model order.
%        .Nsurro: Number of surrogates. Default = 0.
%        .Fs: Sampling frequency
%        .pval
%        .mode: 1-Time; 2-Freq; 3-Both (Default)
%
% To save the data:
%
%   doplot   - true/false (Def: true)
%   dosave   - true/false (Def: false)
%   outpath  - string with path 
%   infosave - string with filename
%
% OUTPUT:
%
% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Feb. 2023; Last revision: 05-Jul-2023

%Change log
% 05-07-2023: Added option to compute only time or only freq GC
% 01-06-2023: Added save options if 2 variables
% 14-02-2023: Plot changed when only 2 variables
% 14-02-2023: Added surrogates and plot them in frequency

%% PARAMETERS

if ~isfield(cnfg,'pval'), cnfg.pval=0.05; end
if ~isfield(cnfg,'dosave'), cnfg.dosave=false; end
if ~isfield(cnfg,'doplot'), cnfg.doplot=true; end
if isfield(cnfg,'outpath')
    if ~strcmp(cnfg.outpath(end),'\')
        cnfg.outpath = [cnfg.outpath '\'];
    end
end
if ~isfield(cnfg,'infosave'), cnfg.infosave=''; end
if ~isfield(cnfg,'mode'), cnfg.mode=3; end

%% PIPELINE

p    = cnfg.p;
w    = cnfg.window;
step = cnfg.step;
Nw   = floor((size(X,2)-w)/step)+1;
Nch  = size(X,1);
fres = cnfg.window;

%% Granger Causality
F = zeros(Nch,Nch,Nw);
wb=waitbar(0,'Computing Granger Causality...');
error_w = [];
for wi = 1:Nw
    Xi = X(:,(wi-1)*step+1:(wi-1)*step+w);
    [A,SIG] = tsdata_to_var(Xi,p);
    [G,info] = var_to_autocov(A,SIG);
    % var_info(info,true)
    if info.error
        error_w = [error_w wi];
    else
        if cnfg.mode==1 || cnfg.mode==3
            F(:,:,wi) = autocov_to_pwcgc(G);
        end
        if cnfg.mode==2 || cnfg.mode==3
            f(:,:,:,wi) = autocov_to_spwcgc(G,fres);
        end
    end
    wb=waitbar(wi/Nw,wb);
end
close(wb)

if cnfg.mode==1 || cnfg.mode==3
    F(:,:,error_w) = [];
    Favg = mean(F,3);
else
    Favg = [];
end

if cnfg.mode==2 || cnfg.mode==3
    if length(error_w)==Nw
        error('Unable to compute GC')
    end
    f(:,:,:,error_w) = [];
    favg = mean(f,4);
else
    favg = [];
end

%% Surrogates

Nsurro = cnfg.Nsurro;

if Nsurro > 0
Fsurro = zeros(Nch,Nch,Nsurro);
wb=waitbar(0,'Computing surrogates...');
error_s = [];
for si = 1:Nsurro
    clear Xsurro
    sp = randperm(size(X,2)-w-1);
    for ch=1:Nch
        Xsurro(ch,:) = X(ch,sp(ch):sp(ch)+w);
    end

    [A,SIG] = tsdata_to_var(Xsurro,p);
    [G,info] = var_to_autocov(A,SIG);
    if info.error
        error_s = [error_s si];
    else
        if cnfg.mode==1 || cnfg.mode==3
            Fsurro(:,:,si) = autocov_to_pwcgc(G);
        end
        if cnfg.mode==2 || cnfg.mode==3
            fsurro(:,:,:,si) = autocov_to_spwcgc(G,fres);
        end
    end 
    wb=waitbar(si/Nsurro,wb);
end
close(wb)

if cnfg.mode==1 || cnfg.mode==3
    Fsurro(:,:,error_s) = [];
    Favg = mean(F,3);
else
    Favg = [];
end

if cnfg.mode==2 || cnfg.mode==3
    fsurro(:,:,:,error_s) = [];
    favg = mean(f,4);
else
    favg = [];
end

Fstats = Favg*0;
fstats = favg*0;
if Nch==2
    Fstats(2,1) = norminv(1-cnfg.pval,mean(squeeze(Fsurro(2,1,:))),std(squeeze(Fsurro(2,1,:))));
    Fstats(1,2) = norminv(1-cnfg.pval,mean(squeeze(Fsurro(1,2,:))),std(squeeze(Fsurro(1,2,:))));

    for fi=1:size(fsurro,3)
        fstats(2,1,fi) = norminv(1-cnfg.pval,mean(squeeze(fsurro(2,1,fi,:))),std(squeeze(fsurro(2,1,fi,:))));
        fstats(1,2,fi) = norminv(1-cnfg.pval,mean(squeeze(fsurro(1,2,fi,:))),std(squeeze(fsurro(1,2,fi,:))));
    end
end
end

%% Plot results
if Nch==2 && cnfg.mode==3
    h=figure; 
    subplot(1,3,1)
    bar([1 2],[Favg(2,1) Favg(1,2)])
    xticks([1 2])
    if isfield(cnfg,'label')
        xticklabels({[cnfg.label{1} '->' cnfg.label{2}],[cnfg.label{2} '->' cnfg.label{1}]})
    else
        xticklabels({'1->2','2->1'})
    end
    if cnfg.Nsurro>0
        hold on,
        plot([0.5 1.5 1.5 2.5],[Fstats(2,1) Fstats(2,1) Fstats(1,2) Fstats(1,2)],'k','LineWidth',2)
    end

    subplot(1,3,2), hold on,
    freq_axis = linspace(0,cnfg.Fs/2,size(favg,3));
    plot(freq_axis,squeeze(favg(2,1,:)),'LineWidth',1)
    ffreq = [freq_axis, fliplr(freq_axis)];
    inBetween = [squeeze(favg(2,1,:))' squeeze(favg(2,1,:))'*0];
    P = fill(ffreq, inBetween, [0, 0.4470, 0.7410],'LineStyle','none');
    alpha(P,0.5)
    axis([0 cnfg.Fs/2 0 max(favg(:)*1.2)])
    if isfield(cnfg,'label')
        title([cnfg.label{1} '->' cnfg.label{2}])
    else
        title('1->2')
    end
    if cnfg.Nsurro>0
        plot(freq_axis,squeeze(fstats(2,1,:)),'k','LineWidth',1)
        inBetween = [squeeze(fstats(2,1,:))' squeeze(fstats(2,1,:))'*0];
        P = fill(ffreq, inBetween, [0.5, 0.5, 0.5],'LineStyle','none');
        alpha(P,0.5)
    end
    xlabel('Frequency (Hz)')

    subplot(1,3,3), hold on,
    plot(freq_axis,squeeze(favg(1,2,:)),'LineWidth',1)
    ffreq = [freq_axis, fliplr(freq_axis)];
    inBetween = [squeeze(favg(1,2,:))' squeeze(favg(1,2,:))'*0];
    P = fill(ffreq, inBetween, [0, 0.4470, 0.7410],'LineStyle','none');
    alpha(P,0.5)
    axis([0 cnfg.Fs/2 0 max(favg(:)*1.2)])
    if isfield(cnfg,'label')
        title([cnfg.label{2} '->' cnfg.label{1}])
    else
        title('2->1')
    end
    if cnfg.Nsurro>0
        plot(freq_axis,squeeze(fstats(1,2,:)),'k','LineWidth',1)
        inBetween = [squeeze(fstats(1,2,:))' squeeze(fstats(1,2,:))'*0];
        P = fill(ffreq, inBetween, [0.5, 0.5, 0.5],'LineStyle','none');
        alpha(P,0.5)
    end
    xlabel('Frequency (Hz)')

    if cnfg.dosave
        if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
        if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
        savefig(h,[cnfg.outpath 'GC' cnfg.infosave])
        saveas(h,[cnfg.outpath 'GC' cnfg.infosave '.png'])
        if cnfg.Nsurro>0
            save([cnfg.outpath 'GC' cnfg.infosave],'Favg','favg','freq_axis','fstats','Fstats')
        else
            save([cnfg.outpath 'GC' cnfg.infosave],'Favg','favg','freq_axis')
        end
    end
    if ~cnfg.doplot, close(h), end

elseif Nch~=2
    if cnfg.mode==1 || cnfg.mode==3
        h=figure; plot_pw(Favg);
        if cnfg.dosave
            if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
            if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
            savefig(h,[cnfg.outpath 'GC_time' cnfg.infosave])
            saveas(h,[cnfg.outpath 'GC_time' cnfg.infosave '.png'])
        end
        if ~cnfg.doplot 
            close(h)
        end
    end
    if cnfg.mode==2 || cnfg.mode==3
        h=figure; plot_spw(favg,cnfg.Fs);
        if cnfg.dosave
            if ~strcmp(cnfg.outpath(end),'\'), cnfg.outpath = [cnfg.outpath '\']; end
            if ~exist(cnfg.outpath,'dir'), mkdir(cnfg.outpath); end
            savefig(h,[cnfg.outpath 'GC_freq' cnfg.infosave])
            saveas(h,[cnfg.outpath 'GC_freq' cnfg.infosave '.png'])
        end
        if ~cnfg.doplot 
            close(h)
        end
    end
    if cnfg.dosave
        save([cnfg.outpath 'GC' cnfg.infosave],'Favg','favg')
    end
end




