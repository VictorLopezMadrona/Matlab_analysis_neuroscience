function data_aw=ftdata2evxic(cnfg,ftdata)

%% Converts data in FT format to a format {Nev,nIC}(ntrials,samples) 
% 
% Syntax:  
%    data_aw=ftdata2evxic(cfg,ftdata);
%
% Inputs:
%   cnfg.eventvalue (Def: 'all')
%
%   ftdata
%
% Outputs:
%   data_aw
%
% See also: evxic2ftdata find_sign_trig aw_signif_comps 

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Oct. 2019; Last revision: 4-Dec-2019

%if ~isfield(cnfg,'eventvalue'), cnfg.eventvalue=unique(ftdata.trialinfo); end
if ~isfield(cnfg,'eventvalue'), cnfg.eventvalue='all'; end
eventvalue=cnfg.eventvalue;

[nch,samples]=size(ftdata.trial{1});

if  strcmp(eventvalue,'all') %do not separate triggers
    data_aw = cell(1,nch);
    for ch=1:nch
        for tr=1:length(ftdata.trial)
            data_aw{ch}(tr,:) = ftdata.trial{tr}(ch,:);
        end
    end
else
    
    %Create matrices with the ICA and trials separately
    data_trlica=cell(1,length(eventvalue));
    for ev=1:length(eventvalue)
        trl_tag=find(ftdata.trialinfo==eventvalue(ev));
        data_trlica{ev}=zeros(nch,samples,length(trl_tag));
        for trl_i = 1:length(trl_tag)
            data_trlica{ev}(:,:,trl_i) = ftdata.trial{trl_tag(trl_i)};
        end
    end

    %Create the struct and analyze the data
    data_aw=cell(length(eventvalue),nch);
    for ev=1:length(eventvalue)
        for nic_i=1:nch
            data_aw{ev,nic_i}=squeeze(data_trlica{ev}(nic_i,:,:))';
        end
    end
end






