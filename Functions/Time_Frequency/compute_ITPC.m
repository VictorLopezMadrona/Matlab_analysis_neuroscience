function itpc=compute_ITPC(F)

%% Compute the ITPC using the output from ft_freqanalysis (freq.fourierspctrm)
% 
% 28/11/2024
% see also ITPC_MEG_SEEG.m

N = size(F,1);
itpc = F./abs(F);     % divide by amplitude
itpc = sum(itpc,1);   % sum angles
itpc = abs(itpc)/N;   % take the absolute value and normalize
itpc = squeeze(itpc); % remove the first singleton dimension

if size(F,2)==1 %Only one channel
    itpc_aux = itpc;
    itpc = [];
    itpc(1,:,:) = itpc_aux;
end
