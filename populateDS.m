% This function populates the struct array with WM load and mean
% voltage data by subject:

function populateDS (filepath, subject, LB, UB)
    global EEG
    Lb = round(-EEG.xmin*EEG.srate+LB*EEG.srate);
    Ub = round(-EEG.xmin*EEG.srate+UB*EEG.srate);
    channel = find(ismember({EEG.chanlocs.labels},'Cz') == 1);

    load([filepath '/' 'TimeCodeEVSdata_' subject '.mat'])
    TextLabels = s.TextLabels;
    clear s
    load('dataframe.mat')
    start = length(s);
    % [txt ',' CWred ',' CWnored,',', SYL,',', CLred,',', CLnored]
    j = length(s);
    for i = 1:length(EEG.epoch)
        tag = str2num(char(EEG.epoch(i).eventtype));
        if length(tag)==6
            j = j + 1;
            s(j).subj = subject;
            s(j).txt_no = tag(1);
            s(j).lang = TextLabels{tag(1), 2};
            s(j).CWred = tag(2);
            s(j).CWnored = tag(3);
            s(j).SYL = tag(4);
            s(j).CLred = tag(5);
            s(j).CLnored = tag(6);
            s(j).mu_stat = mean(EEG.data(channel, Lb:Ub, i));
            s(j).erp = EEG.data(channel, :, i);
            s(j).epochN = i;
            s(j).fullEp = squeeze(EEG.data(:, :, i));
        end
    end
save('dataframe.mat', 's')