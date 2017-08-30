% This script chunks continuous EEG data into epochs of fixed size
% time-locked to probe onset times. See the manuscript for details.

% clean up junk events:
index = find(strcmp({EEG.event.type}, 'boundary') == 1);
EEG.event(index) = [];
index = find(strcmp({EEG.event.type}, 'S 10') == 1);
EEG.event(index) = [];
index = find(strcmp({EEG.event.type}, 'S 11') == 1);
EEG.event(index) = [];
index = find(strcmp({EEG.event.type}, 'S  0') == 1);
EEG.event(index) = [];

% create 600-ms epochs and save the new dataset:
% create epochs with 100 ms pre-stimulus onset and 500 ms post-stimulus:
EEG = pop_epoch( EEG, unique({EEG.event.type}), [-0.1 0.5], 'epoched', 'yes');

% check the dataset:
EEG = eeg_checkset( EEG );

% remove baseline (based on a 100-ms pre-stimulus period):
EEG = pop_rmbase( EEG, [-100 0]);

% check the dataset:
EEG = eeg_checkset( EEG );

% create a new EEGlab dataset:
name = EEG.setname;
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', name);

% check again:
EEG = eeg_checkset( EEG );

% refresh:
eeglab redraw;