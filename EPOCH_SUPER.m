index = find(strcmp({EEG.event.type}, 'boundary') == 1);
EEG.event(index) = [];

index = find(strcmp({EEG.event.type}, 'S 10') == 1);
EEG.event(index) = [];

index = find(strcmp({EEG.event.type}, 'S 11') == 1);
EEG.event(index) = [];

index = find(strcmp({EEG.event.type}, 'S  0') == 1);
EEG.event(index) = [];

EEG = pop_epoch( EEG, unique({EEG.event.type}), [-0.1 0.5], 'epoched', 'yes');
EEG = eeg_checkset( EEG );
EEG = pop_rmbase( EEG, [-100 0]);
EEG = eeg_checkset( EEG );
name = EEG.setname;
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', name);
EEG = eeg_checkset( EEG );
eeglab redraw;