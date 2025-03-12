
function [ftdata, elec_pick] = VBMEG2FT(filename)

% Load a VBMEG file and convert it to FieldTrip. Continuous data, no trials 
%
% Victor J. Lopez-Madrona
% 10/03/2025

vbdata = load(filename);

ftdata = [];
%ftdata.hdr
ftdata.label = vbdata.MEGinfo.ChannelInfo.Name;
ftdata.trial{1} = vbdata.bexp;
ftdata.fsample = vbdata.MEGinfo.SampleFreq;
ftdata.time{1} = linspace(0,size(ftdata.trial{1},2)/ftdata.fsample,size(ftdata.trial{1},2));
elec_pick = vbdata.pick;




