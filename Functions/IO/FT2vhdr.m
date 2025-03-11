function FT2vhdr(ftdata,filepath_vhdr)

% Convert FieldTrip data into .vhdr and save it.
%
% Victor J. Lopez-Madrona
% 10/03/2025

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% WORK IN PROGRESS %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% So far only works with continuous data from FieldTrip file without
% markers

[pathname, filename] = fileparts(filepath_vhdr);
data = ftdata.trial{1};
labels = ftdata.label;
marker = [];
Fsample = ftdata.fsample;

dyn_write_brainvision(pathname, filename, data, labels, marker, Fsample);