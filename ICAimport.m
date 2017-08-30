% The script loads previously computed ICA information and pastes it into
% the current dataset.

filepath = ['/Users/RomanKoshkin/Documents/MATLAB/EEG/' EEG.subject '_EEG/'];
load([filepath 'ICA.mat']) % load ICA weights

EEG.icaact = icaact;
EEG.icawinv = icawinv;
EEG.icaweights = icaweigts;
EEG.icachansind = icachansind;
EEG.icasphere = icasphere;
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw

% open the dialog for reviewing the independent components, their topographies
% and marking them for rejection:
pop_selectcomps(EEG,[1:length(EEG.icaweights)])