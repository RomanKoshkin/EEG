
%% PARAMETERS:

% ASR paramters:
arg_flatline = 5;
arg_highpass = [0.25 0.75];
arg_channel = 0.8;
arg_noisy = 'off';
arg_burst = 3.5;
arg_window = 0.3;
    
% vector of estimators:
interpolate = 1; % 1 use linear interpolation to estimate WM load at probe onsets
                 % 0 use on interpolation, assume WM load step changes at
                 % word offsets

% defines the boundary values. Choose 'quant' for quantiles, 'mean' or 'median'
cutoff = 'quant';

% select two of the following: %[txt ',' CWred ',' CWnored,',', SYL,',', CLred,',', CLnored]
% by specifying their ordinal number:
est_vec = [5 2];

global EST;

% tweak condition boundaries in the cell 'EST'. These condition boundaries
% only apply if you are interested in imposing fixed cutoffs between high
% and low WM load regardless of the subject or direction of SI. See the
% manuscript for details:
load('params.mat');

global EEG;
subjects =  {'KOK', 'GRU', 'ROM', 'SHE' 'KOS', 'ELT', 'POG', 'BUL', 'KOZ'};

% genetic algorithm for EEG cleaning (not used). Set true,
% if you want to try it. It is just a proof of concept and we did not use
% it in the paper.
genetic = true;

% Lower and upper bounds (in seconds, post stimulus),for calculating mean 
% votage in the N1 interval:
LB = 0.12;
UB = 0.18;

% make an empty data structure:
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

% this part of the code generates a data struct needed for the legacy
% grand_avg_SUPER.m to plot grand average ERPs.Use NEW_GRAND_AVERAGE.m 
% instead.
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

% start the timer:
timerx = clock;
tic

% start EEGlab:
eeglab

for cycle = 1:length(subjects)
    
    Setname = [subjects{cycle} '_RAWset3_rec'];
    OUTfilename = [Setname '_ASR.set'];
    INfilename =  [subjects{cycle} '_RAWset3.set'];
    filepath = ['/Users/RomanKoshkin/Documents/MATLAB/EEG/' subjects{cycle} '_EEG'];
    
 if exist([filepath '/' OUTfilename], 'file') == 2
     % if the file cleaned with ASR already exists, skip the ASR stage:
    disp(['file ' OUTfilename ' already exists in the path ' filepath '. No need to re-create'])
     EEG = pop_loadset([filepath '/' OUTfilename]);
     eeglab redraw
 
 else
     % clean with ASR:
    EEG = pop_loadset([filepath '/' INfilename]);
    eeglab redraw
    EEG = clean_rawdata(EEG, arg_flatline, arg_highpass, arg_channel, arg_noisy, arg_burst, arg_window);
    EEG = eeg_checkset( EEG );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', Setname);
    eeglab redraw
    EEG = eeg_checkset( EEG );
    EEG = pop_saveset(EEG,'filename', OUTfilename, 'filepath', filepath, 'check', 'on', 'savemode', 'onefile');
 end
    if length([EEG.chanlocs]) < 30 
        % if the number of channels is less than 30, interpolate the
        % missing ones:
        interp_chan
    end
    
    % Change event codes to match WM load estimates:
    EVS_corr_SUPER(interpolate)
    
    % Epoch the current subject's dataset:
    EPOCH_SUPER
    
    % populateDS(filepath, subjects{cycle}, LB, UB)
    
    % plot ERPs for each of the subjects:
    ERPimg_SUPER(filepath, subjects{cycle}, est_vec, genetic, cutoff, LB, UB)
    
    ALLCOM = []; ALLEEG = []; ALLERP = []; CURRENTSET = 0; CURRENTSTUDY = 0; EEG = [];
    eeglab redraw
    disp (['end of cycle' subjects{cycle}])
    
end

% stop the timer
toc

% legacy function to plot grand average ERPs. Use NEW_GRAND_AVERAGE.m instead.
grand_avg_SUPER(subjects, est_vec, LB, UB)

% report time elapsed:
clock - timerx