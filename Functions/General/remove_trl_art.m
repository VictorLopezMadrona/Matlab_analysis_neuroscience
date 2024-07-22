
function [ftdata]=remove_trl_art(cnfg,ftdata)

%% Removes trials with artefacts in MEG data
% 
% Syntax:  
%    [ftdata]=remove_trl_art(cfg,ftdata);
%
% Inputs:
%   cfg.filename     = Path with the .mrk file.
%      .marker_name  = [Opt] Name of the markers with artefacts. Def='artefact'
%      .marker_value = [Opt] Value of the markers with artefacts. Def=-99
%        
%   ftdata = struct as the one given by ft_preprocessing.
%
% Outputs:
%   ftdata = same as the original but without trials with artefacts
%   ftdata.trialsremoved = Trials that were removed in the original data.    
%
% See also: loadMEGdata

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Nov. 2019; Last revision: 27-Nov-2019

%% Default values

%Marker by name. 
if isfield(cnfg,'marker_name') 
    if iscell(cnfg.marker_name)
        name=cnfg.marker_name;
    else
        name{1}=cnfg.marker_name;
    end
else 
    name={'artefact'}; 
end

%Marker by value
if isfield(cnfg,'marker_value') 
    if iscell(cnfg.marker_value)
        name=cnfg.marker_value;
    else
        name{1}=cnfg.marker_value;
    end
else 
    value={'-99'}; 
end

%Check if the file with markers exist
if isfield(cnfg,'filename'), mark_fullfile=cnfg.filename;
else, error('The name of the marker file must be specified'), end
if ~exist(mark_fullfile,'file')
    warning('No file with markers found. Trials with artefacts will not be removed')
    return
end

%% Workflow

fileID = fopen(mark_fullfile);
MARKERS = textscan(fileID,'%s %s %s %s');
fclose(fileID);
for i=1:4
    MARKERS{i}=MARKERS{i}(2:end); %Remove the first element
end

% disp('Markers found in the file (by name):') 
% disp(unique(MARKERS{1})) %Markers in the file
% disp('Markers found in the file (by value):') 
% disp(unique(MARKERS{2})) %Markers in the file

p_art = zeros(length(MARKERS{1}),length(name)+length(value));
for n=1:length(name)
    p_art(:,n) = ismember(MARKERS{1},name{n}); %find artefacts
end
for n=1:length(value)
    p_art(:,n+length(name)) = ismember(MARKERS{2},num2str(value{n})); %find artefacts
end
p_art = sum(p_art,2)>0;

if sum(p_art)==0
    disp('No artefacts found in the file')
else
    t_art = [str2double(MARKERS{3}(p_art)) str2double(MARKERS{3}(p_art))+str2double(MARKERS{4}(p_art))];
    sample_art = t_art*ftdata.fsample;
    rmv_trl = zeros(1,length(ftdata.trial));
    for nt=1:length(ftdata.trial)
        trl_art(:,1) = ftdata.sampleinfo(nt,1)>sample_art(:,1); %artefacts start before ini trial 
        trl_art(:,2) = ftdata.sampleinfo(nt,1)<sample_art(:,2); %artefacts end after ini trial 
        trl_art(:,3) = ftdata.sampleinfo(nt,1)<sample_art(:,1); %artefacts start after ini trial 
        trl_art(:,4) = ftdata.sampleinfo(nt,2)>sample_art(:,1); %artefacts start before end trial 
        %Test first two conditions or last two. Both imply an overlap
        %between trial and artefact
        overlp = (trl_art(:,1) & trl_art(:,2)) | (trl_art(:,3) & trl_art(:,4));
        if sum(overlp)>0, rmv_trl(nt)=1; end
    end
    
    if sum(rmv_trl)==0
        disp('No artefacts found during trials')
    else
        ftdata.time = ftdata.time(rmv_trl==0);
        ftdata.trial = ftdata.trial(rmv_trl==0);
        ftdata.sampleinfo = ftdata.sampleinfo(rmv_trl==0,:);
        ftdata.trialinfo = ftdata.trialinfo(rmv_trl==0);

        warning([num2str(sum(rmv_trl)) ' trials with artefacts were removed'])
    end
    [~,ftdata.trialsremoved] = find(rmv_trl==1);
end
        

        
        
    
    