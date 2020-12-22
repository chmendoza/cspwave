function EEG = mat2EEGstruct(data, srate)
% data is CxN, C: # of channels. N: # of samples. srate is sampling rate in
% samples/s.

%% Initialize EEG struct (EEGLAB and PREP pipeline)
EEG = eeg_emptyset;
EEG.srate = srate;
%% Populate EEG struct
EEG.data = data;
EEG.nbchan = size(EEG.data, 1);
EEG.pnts = size(EEG.data, 2);
t = linspace(0, EEG.pnts/srate, EEG.pnts);
EEG.times = t;

end