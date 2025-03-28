function ftdata=interpol_artefact(cnfg,ftdata)

%% Function to do a basic interpolation in a time-window with artefact
%
% 19/03/2025

Ntr = length(ftdata.trialinfo);
Nch = length(ftdata.label);
tini = find(ftdata.time{1}>=cnfg.time(1),1);
tend = find(ftdata.time{1}>=cnfg.time(2),1);

for tri=1:Ntr
    for chi=1:Nch
        dur = tend-tini+1;
        noise_trial = ftdata.trial{tri}(chi,tend)-ftdata.trial{tri}(chi,tini);
        new_segment = linspace(ftdata.trial{tri}(chi,tini),ftdata.trial{tri}(chi,tend),dur) + (randn(1,dur)*noise_trial/5);
        ftdata.trial{tri}(chi,tini:tend) = new_segment;
    end
end

    
    

