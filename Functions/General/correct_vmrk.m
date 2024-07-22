function correct_vmrk(filename,marker_name)

%% For vmrk files, the responses should be defined as 'R' and stimulus as 'S'
% This function will replace those elements in the .vmrk file

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Sep. 2021; Last revision: 20-Sep-2021

if nargin == 1
    marker_name{1,1} = 'MARQUEUR'; marker_name{1,2} = 'S';
    marker_name{2,1} = 'Stimulus'; marker_name{2,2} = 'S';
    marker_name{3,1} = 'RESPONSE'; marker_name{3,2} = 'R';
end

fileID = fopen(filename);
lines_mrk = textscan(fileID,'%s','delimiter','\n');
lines_mrk = lines_mrk{1};
fclose(fileID);

for i=1:size(marker_name,1)
    for l=6:length(lines_mrk)
        lines_mrk{l} = strrep(lines_mrk{l},marker_name{i,1},marker_name{i,2});
    end
end

% Write results in a new marker file
fid=fopen(filename,'w');
for i=1:length(lines_mrk)
    fprintf(fid,[lines_mrk{i} '\n']);
end
fclose(fid);


