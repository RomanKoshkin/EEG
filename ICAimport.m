filepath = ['/Users/RomanKoshkin/Documents/MATLAB/EEG/' EEG.subject '_EEG/'];

load([filepath 'ICA.mat']) % load ICA weights


EEG.icaact = icaact;
EEG.icawinv = icawinv;
EEG.icaweights = icaweigts;
EEG.icachansind = icachansind;
EEG.icasphere = icasphere;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw
pop_selectcomps(EEG,[1:length(EEG.icaweights)])