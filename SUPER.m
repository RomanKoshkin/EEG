    % ASR paramters:
arg_flatline = 5;
arg_highpass = [0.25 0.75];
arg_channel = 0.8;
arg_noisy = 'off';
arg_burst = 3.5;
arg_window = 0.3;
    % END ASR PARAMETERS =+++++++++++++++++++++++++++++++++++++++++++++++
    
% vector of estimators:
interpolate = 1;
cutoff = 'quant'; % 'mean' or 'median'
est_vec = [5 2];  %[txt ',' CWred ',' CWnored,',', SYL,',', CLred,',', CLnored]
global EST; load('params.mat') % tweak condition boundaries here
global EEG;
subjects =  {'KOK', 'GRU', 'ROM', 'SHE' 'KOS', 'ELT', 'POG', 'BUL', 'KOZ'};
genetic = false;
LB = 0.12; % seconds. Lower bound for measuring the mean votage in the N1 interval
UB = 0.18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear the general dataframe:
h = struct(...
    'subj', {'a', 'a'},...
    'txt_no', {'a', 'a'},...
    'lang', {'a', 'a'},...
    'CWred', 0,...
    'CWnored', 0,...
    'SYL', 0,...
    'CLred', 0,...
    'CLnored', 0,...
    'mu_stat', 0,...
    'erp', 0);
save('dataframe.mat', 'h')
clear h
%clear all the fields in AV:
load('AV.mat')
I = fields(AV);
J = fields(AV.SHE);
K = fields(AV.SHE.both);
L = fields(AV.SHE.both.CW);
M = fields(AV.SHE.both.CW.C1);
for i = 1:length(I)
    for j = 1:length(J)
        for k = 1:length(K)
            for l = 1:length(L)
                for m = 1:length(M)
                    AV.(I{i}).(J{j}).(K{k}).(L{l}).(M{m}) = 0;
                end
            end
        end
    end
end
save('AV.mat', 'AV')
clear I J K L M i j k l m AV
% END clear AV

timerx = clock;
tic
eeglab
for cycle = 1:length(subjects)
    
%     Setname = [subjects{cycle} '_RAWset3_ASR('...
%         num2str(arg_channel) ',' arg_noisy ',' ...
%         num2str(arg_burst) ',' num2str(arg_window) ')'];
    Setname = [subjects{cycle} '_RAWset3_rec'];
    OUTfilename = [Setname '_ASR.set'];
    INfilename =  [subjects{cycle} '_RAWset3.set'];
    
    filepath = ['/Users/RomanKoshkin/Documents/MATLAB/EEG/' subjects{cycle} '_EEG'];
 if exist([filepath '/' OUTfilename], 'file') == 2
    disp(['file ' OUTfilename ' already exists in the path ' filepath '. No need to re-create'])
     EEG = pop_loadset([filepath '/' OUTfilename]);
     eeglab redraw
 else
    EEG = pop_loadset([filepath '/' INfilename]);
    eeglab redraw
    EEG = clean_rawdata(EEG, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst, arg_window);
    EEG = eeg_checkset( EEG );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', Setname);
    eeglab redraw
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename', OUTfilename, 'filepath', filepath, 'check', 'on', 'savemode', 'onefile');
 end
    if length([EEG.chanlocs]) < 30 % if the number of channels is less than 30, interpolate missing
        interp_chan
    end
    EVS_corr_SUPER(interpolate) % 1 - Estimate interpolated load values at probe onset, 0 - do not interpolate
    EPOCH_SUPER
    % populateDS(filepath, subjects{cycle}, LB, UB)
    ERPimg_SUPER(filepath, subjects{cycle}, est_vec, genetic, cutoff, LB, UB)
    
    ALLCOM = []; ALLEEG = []; ALLERP = []; CURRENTSET = 0; CURRENTSTUDY = 0; EEG = [];
    eeglab redraw
    disp (['end of cycle' subjects{cycle}])
    
end
toc
grand_avg_SUPER(subjects, est_vec, LB, UB)
clock - timerx