function dyn_write_brainvision(pathname, filename, data, labels, marker, Fsample)

% Function from Dynamap Team, INS, Marseille, France
%
% TO-DO
% Check the format of the data: EEG, MEG...

nbr_mark = length(marker);

%Init value
Nchan = length(labels);
Tsample = (1/Fsample)*10^6;
Bformat = 'IEEE_FLOAT_32';

%Create the filename for the chunks
vhdr_file = fullfile(pathname, [filename '.vhdr']);
vmrk_file = fullfile(pathname, [filename '.vmrk']);
eeg_file = fullfile(pathname, [filename, '.eeg']);

%get the data for the binary file
if isstruct(data)
    data_binary = zeros(Nchan, length(data(1).data));
    for n=1:Nchan
        data_binary(n, :) = data(n).data(:);
    end
else
    data_binary = data;
end
fileID = fopen(eeg_file, 'w');
fwrite(fileID, data_binary, 'float32', 'ieee-le');
fclose(fileID);

%Create the vhdr and the vmrk file
fHeader = fopen(vhdr_file, 'w+t');
fMarker = fopen(vmrk_file, 'w+t');
%Create the string value
Header = 'Brain Vision Data Exchange Header File Version 1.0\n; File created by AnyWave\n\n[Common Infos]\n';
MHeader = 'Brain Vision Data Exchange Marker File Version 1.0\n; File created by AnyWave\n\n[Common Infos]\n';

%Create the header file .vhdr
fprintf(fHeader, Header);
fprintf(fHeader, 'DataFile=%s\n', strcat(filename, '.eeg'));
fprintf(fHeader, 'MarkerFile=%s\n', strcat(filename, '.vmrk'));
fprintf(fHeader, 'DataFormat=BINARY\n');
fprintf(fHeader, 'DataType=TIMEDOMAIN\n');
fprintf(fHeader, 'DataOrientation=MULTIPLEXED\n');
fprintf(fHeader, 'NumberOfChannels=%d\n', Nchan);
fprintf(fHeader, 'SamplingInterval=%.3f\n', Tsample);

fprintf(fHeader, '\n[Binary Infos]\n');
fprintf(fHeader, 'BinaryFormat=%s\n', Bformat);

fprintf(fHeader, '\n[Channel Infos]\n');
for i=1:Nchan
    %fprintf(fHeader, 'Ch%d=%s,,1,µV\n', i, labels{i});
    fprintf(fHeader, 'Ch%d=%s,,1,T\n', i, labels{i});
end

%Create the marker file
fprintf(fMarker, MHeader);
fprintf(fMarker, 'DataFile=%s\n', strcat(filename, '.eeg'));
fprintf(fMarker, '\n[Marker Infos]\n');
if isstruct(marker)
    fields = fieldnames(marker);
    if any(strcmp(fields, 'label')) && any(strcmp(fields, 'position')) && any(strcmp(fields, 'duration'))
        for i=1:nbr_mark
            fprintf(fMarker, 'Mk%d=Comment,%s,%d,%d,0\n',i, marker(i).label, floor(marker(i).position*Fsample), marker(i).duration*Fsample);
        end
    elseif any(strcmp(fields, 'value')) && any(strcmp(fields, 'sample')) && any(strcmp(fields, 'duration'))
        for i=1:nbr_mark
            fprintf(fMarker, 'Mk%d=Comment,%s,%d,%d,0\n',i, marker(i).value, marker(i).sample, marker(i).duration);
        end
    end
end
fclose(fMarker);
fclose(fHeader);
end
