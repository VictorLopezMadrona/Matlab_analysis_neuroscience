function src_raw = conv_source_to_raw(source)
%% Convert Field Trip's source analysis structure to raw type of data
%
% Syntax:  
%    src_raw = conv_source_to_raw(source)
%
% Inputs:
%    source - Source analysis type of structure, the output of ft_sourceanalysis
%
% Outputs:
%    src_raw - The converted source structure to am equivalent raw type of
%              data, see ft_preprocessing
%
% See also: ft_sourceanalysis,  ft_preprocessing

% Author: Christos Papageorgakis <christos.papageorgakis@gmail.com>
% License: BSD (3-clause)
% Sept. 2019; Last revision: 06-Sept-2019

% memopt = true;
% if memopt
%     avg_names = fieldnames(source.avg);
%     % remove mom
%     avg_names(find(ismember(avg_names, 'mom'))) = [];
%     for idx=1:length(avg_names)
%         source.avg = rmfield(source.avg, avg_names{idx});
%     end
% end

%% Convert source to raw type of field trip data
src_raw = [];
src_raw.time = {source.time};
src_raw.fsample = length(source.time)/source.time(end);

%%% cell2mat is slower and demands more memory than the following for loop
% src_raw.trial = {cell2mat(source.avg.mom(source.inside))};

%%% Convert row cells to matrix
nbr_samples = length(source.time);
nbr_channels = nnz(source.inside);
src_raw.trial = nan(nbr_channels, nbr_samples);
inside_idcs = find(source.inside);
for curr_ch = 1:nbr_channels
    src_raw.trial(curr_ch, :) = source.avg.mom{inside_idcs(curr_ch)};
    source.avg.mom{inside_idcs(curr_ch)} = [];
end
src_raw.trial = {src_raw.trial};

src_raw.sampleinfo = [1, length(src_raw.trial{1})];
src_raw.label =  strsplit(num2str(1:length(source.inside)))';
src_raw.label =  src_raw.label(source.inside);
% cnfg.labels = strcat('Ch_', cnfg.labels);

end


