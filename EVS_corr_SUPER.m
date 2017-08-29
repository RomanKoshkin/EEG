function EVS_corr_super(interpolate)
head = {'CWnored', 'onset', 'offset', 'SYL', 'CWred',...
    'CLnored', 'CLred', 'probe_time', 'CWnored_interp',...
    'SYL_interp', 'CWred_interp', 'CLnored_interp', 'CLred_interp'};

global EEG ALLEEG CURRENTSET
subject = EEG.subject;
path = ['/Users/RomanKoshkin/Documents/MATLAB/EEG/', subject, '_EEG/'];
load([path 'TimeCodeEVSdata_' subject '.mat'])

if isfield(s, 'TextLabels') == 1 % remove the unnecessary field 'TextLabels' if it exists
    s = rmfield (s, 'TextLabels');
end

for text_no = 1:length(fieldnames(s))
        a = 0;
        txt = num2str(text_no);
        for i = 1:length(EEG.event)
            if EEG.event(i).type == txt
               a(i) = EEG.event(i).latency/EEG.srate;
            else
               a(i) = NaN;
            end
        end

    tempmt = s.(['text' txt]);
    edges = [tempmt(:,2); tempmt(end,3)];
    data = a;
    Y = discretize(data,edges);
    x = [a;Y; NaN(2, length(Y))];
    for i = 1:length(edges)-1
        x(3,find(x(2,:)==i)) = edges(i);
        x(4,find(x(2,:)==i)) = edges(i+1);
    end
    %%%%%%%%%%%%%%%%%%%%% now interpolate shit: %%%%%%%%%%%%%%%%
    if interpolate == 1
        for i = 1:length(x)
            C = x(1,i);
            timeONvec = tempmt(:,2);
            timeOFFvec = tempmt(:, 3);
            addr = find (timeONvec <= C & C < timeOFFvec);
            if addr >= size(tempmt,1)
                break
            end

            A = timeONvec(addr);
            B = timeOFFvec(addr);

            tempmt(addr,8) = C;

            % CWnored_interp
            tempmt(addr,9) = (C-A)/(B-A)*(tempmt(addr+1, 1) - tempmt(addr, 1)) + tempmt(addr, 1);
            % SYL_interp
            tempmt(addr,10) = (C-A)/(B-A)*(tempmt(addr+1, 4) - tempmt(addr, 4)) + tempmt(addr, 4);
            % CWred_interp
            tempmt(addr,11) = (C-A)/(B-A)*(tempmt(addr+1, 5) - tempmt(addr, 5)) + tempmt(addr, 5);
            % CLnored_interp
            tempmt(addr,12) = (C-A)/(B-A)*(tempmt(addr+1, 6) - tempmt(addr, 6)) + tempmt(addr, 6);
            % CLnored_interp
            tempmt(addr,13) = (C-A)/(B-A)*(tempmt(addr+1, 7) - tempmt(addr, 7)) + tempmt(addr, 7);
        end

        % T = array2table(tempmt, 'VariableNames', head);
        % plotinterp(tempmt, subject, text_no)

        for i = 1:length(data)
            idx = x(2,i);
            if ~isnan(idx)
                CWnored = num2str(tempmt(idx, 9));
                SYL = num2str(tempmt(idx, 10));
                CWred = num2str(tempmt(idx, 11));
                CLnored = num2str(tempmt(idx, 12));
                CLred = num2str(tempmt(idx, 13));

                EEG.event(i).type = [txt ',' CWred ',' CWnored,',', SYL,',', CLred,',', CLnored]; % '@' num2str(s.text1(idx, 2)) '-' num2str(s.text1(idx, 3))];

            end
        end
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:length(data)
            idx = x(2,i);
            if ~isnan(idx)
                CWred = num2str(tempmt(idx, 1));
                CWnored = num2str(tempmt(idx, 5));
                SYL = num2str(tempmt(idx, 4));
                CLred = num2str(tempmt(idx, 7));
                CLnored = num2str(tempmt(idx, 6));

                % modify the event triggers to contain the five different estimators:
                EEG.event(i).type = [txt ',' CWred ',' CWnored,',', SYL,',', CLred,',', CLnored]; % '@' num2str(s.text1(idx, 2)) '-' num2str(s.text1(idx, 3))];

            end
        end
    end
end
disp (['number of texts in the EVS file: ' num2str(text_no) ' Subject: ' subject])

qq = [EEG.subject,'_RAWset4'];
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', qq);
EEG = eeg_checkset( EEG );
eeglab redraw
end