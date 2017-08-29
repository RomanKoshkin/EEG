typerange = {'1','2','3','4', '5', '6', '7', '8'};
subjects =  {'GRU', 'ROM', 'SHE' 'KOS', 'ELT', 'POG', 'BUL', 'KOZ', 'KOK'};
for cycle = 1:length(subjects)        
    Setname = [subjects{cycle} '_RAWset3'];
    OUTfilename = [Setname, '_rec', '.set'];
    INfilename =  [subjects{cycle} '_RAWset3.set'];
    filepath = ['/Users/RomanKoshkin/Documents/MATLAB/EEG/' subjects{cycle} '_EEG'];

    EEG = pop_loadset([filepath '/' INfilename]);
    EEG = pop_eegfiltnew(EEG, 0.25, 30);
    EEG = eeg_checkset( EEG );
    name = EEG.setname;
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', name);
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename', OUTfilename, 'filepath', filepath, 'check', 'on', 'savemode', 'onefile');
    
    OUTEEG = pop_rmdat( EEG, typerange , [-0.5 0.5], 0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, OUTEEG, CURRENTSET, 'setname', 'temp');
    
%     EEG = pop_runica(EEG, 'icatype', 'binica', 'extended', 1);
%     [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%     icaact = EEG.icaact;
%     icawinv = EEG.icawinv;
%     icaweigts = EEG.icaweights;
%     icachansind = EEG.icachansind;
%     icasphere = EEG.icasphere;
%     % save('ICA.mat', 'icaact',  'icawinv', 'icaweigts', 'icachansind', 'icasphere')
%     ALLEEG = pop_delset(ALLEEG, 2);
%     eeglab redraw
%     
%     EEG.icaact = icaact;
%     EEG.icawinv = icawinv;
%     EEG.icaweights = icaweigts;
%     EEG.icachansind = icachansind;
%     EEG.icasphere = icasphere;
%     [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%     
%     EEG = pop_saveset(EEG,'filename', OUTfilename, 'filepath', filepath, 'check', 'on', 'savemode', 'onefile');
%     clear OUTEEG icaact icawinv icaweigts icachansind icasphere
%     ALLEEG = pop_delset(ALLEEG, 1);
end

eeglab redraw