% this function interpolates missing channels. In some subjects the signal
% was very noisy on temporal channels FT9 FT10 T7 T8 and had to be removed.

function interp_chan()
    global EEG ALLEEG CURRENTSET
    load('chanlocs.mat');
    filepath = EEG.filepath;
    filename = EEG.filename;
    EEG = pop_interp(EEG, chanlocs, 'spherical');
    name = [EEG.setname '_interp'];
    EEG = eeg_checkset( EEG );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', name);
    eeglab redraw
    EEG = pop_saveset(EEG,'filename', filename, 'filepath', filepath, 'check', 'on', 'savemode', 'onefile');
end