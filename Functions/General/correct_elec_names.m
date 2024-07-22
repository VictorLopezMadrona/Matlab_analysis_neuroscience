function label=correct_elec_names(label,flag_bip)

if nargin==1
    flag_bip=[];
end

%Add '_' to separate name and number
%Eliminate multiple '_' in SEEG names
%Change ' by p
%If bipolar, remove reference.

if strcmp(flag_bip,'bip')
for l=1:length(label)
    p_refbip=strfind(label{l},'-');
    label{l}(p_refbip:end)=[];
end
end

for l=1:length(label)
    isn = str2double(label{l}(end-1:end));
    if ~isnan(isn)
        label{l} = [label{l}(1:end-2) '_' label{l}(end-1:end)];
    else
        if ~isnan(str2double(label{l}(end)))
            label{l} = [label{l}(1:end-1) '_' label{l}(end)];
        else
            warning('Name and number not identified in electrode name')
        end
    end
end

for l=1:length(label)
    flag=1;
    while flag==1
        label_old = label{l};
        label{l}=regexprep(label{l},'__','_');
        if strcmp(label{l},label_old)
            flag=0;
        end
    end
end

for l=1:length(label)
   label{l}=regexprep(label{l},"'",'p');
end 
    
    
    
    