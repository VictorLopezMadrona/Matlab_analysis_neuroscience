
function [ftdata, elec_pick] = VBMEG2FT(filename)

%% Load a VBMEG file and convert it to FieldTrip. Continuous data, no trials 
%

% Author: Victor Lopez Madrona <v.lopez.madrona@gmail.com>
% License: BSD (3-clause)
% Mar. 2025; Last revision: 10-Apr-2025


warning('''elec_pick'' output will be removed in future versions. Use ftdata.grad.coilpos instead')

vbdata = load(filename);

ftdata = [];
%ftdata.hdr
ftdata.label = vbdata.MEGinfo.ChannelInfo.Name;
ftdata.trial{1} = vbdata.bexp;
ftdata.fsample = vbdata.MEGinfo.SampleFreq;
ftdata.time{1} = linspace(0,size(ftdata.trial{1},2)/ftdata.fsample,size(ftdata.trial{1},2));

grad=[];
grad.chanpos = vbdata.pick;
grad.chanori = vbdata.Qpick;
grad.chantype = repmat({vbdata.MEGinfo.Measurement}, vbdata.MEGinfo.Nchannel, 1);
grad.chanunit = repmat({'T'}, vbdata.MEGinfo.Nchannel, 1);
grad.coilpos = vbdata.pick;
grad.coilori = vbdata.Qpick;
grad.unit = 'm';
%grad.coordsys = 'acpc';
grad.coordsys = vbdata.CoordType;
grad.label = vbdata.MEGinfo.MEGch_name;
ftdata.grad = grad;

elec_pick = vbdata.pick;


