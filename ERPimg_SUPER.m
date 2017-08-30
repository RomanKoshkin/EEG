% This legacy function plots ERPs for a specified subject.
% Unlike the NEW_GRAND_AVERAGE.m function, it does not compute
% subject- and direction-specific cutoffs. The cutoff values for every WM
% load estimator are preset in the params.mat file.

% Parameters:
% filepath      home directory of the current subject.
% subject       subject code ('KOS', 'BUL', ..., etc.)
% est_vec       a 2-vector of estimators. For example: [2 3] tells the function
%               that we want to compare CWred and CLred ([txt CWred CWnored SYL CLred
%               CLnored])
% genetic       set this flag to 'true' if you want to use the function
%               GA_SUER.m to filter out bad epochs (EXPERIMENTAL).
%               Otherwise set it to zero.
% data sep      criterion for splitting the epochs based on different
%               levels of WM load. 
%                   quant - split that data along the quantiles specified
%                   by LB (boundary between low and medium WM load) and UB
%                   (quantile between medium and high WM load).
% LB and UB     see above.

function ERPimg_SUPER(filepath, subject, est_vec, genetic, data_sep, LB, UB)
    global EEG EST
    
    % electrode of interest:
    channel = 'Cz';
    channel = find(ismember({EEG.chanlocs.labels}, channel) == 1);
    
    % load WM load data:
    load([filepath '/' 'TimeCodeEVSdata_' subject '.mat'])
    X = {EEG.event.type}';
    for i = 1:length(X)
         temp = str2num(X{i});
         if length(temp)==6
             X1(i,:) = [temp EEG.event(i).epoch];
         else
             X1(i,:) = [NaN NaN NaN NaN NaN NaN NaN];
         end
    end

    for Dspl = [0 2 4]
        switch Dspl
            case 0
                txtNo = [1 2 3 4 5 6 7 8];        % text labels
                language = 'both';
                counter = 0;
            case 2
                txtNo = find(strcmp(s.TextLabels(:,2),'er'));
                language = 'ER';
            case 4
                txtNo = find(strcmp(s.TextLabels(:,2),'RE'));
                language = 're';
        end
        
        % get epochs from the right texts (EN, RU, or both):
        mask = arrayfun(@(x) find(X1(:,1)==x), txtNo, 'UniformOutput', false);
        mask = cat(1, mask{:}); % vectorize the cell array
        masked = X1;
        masked(setdiff(1:length(masked), mask), :) = NaN;

        %[txt CWred CWnored SYL CLred CLnored]
        if counter == 0
            load('dataframe.mat')
            start = length(h);
            for jk = 1:length(masked)
                if isnan(masked(jk,1)) == 0
                    counter = counter + 1;

                    h(start+counter).subj = subject;
                    h(start+counter).txt_no = masked(jk,1);
                    h(start+counter).lang = s.TextLabels{masked(jk,1), 2};
                    h(start+counter).CWred = masked(jk,2);
                    h(start+counter).CWnored = masked(jk,3);
                    h(start+counter).SYL = masked(jk,4);
                    h(start+counter).CLred = masked(jk,5);
                    h(start+counter).CLnored = masked(jk,6);
                    h(start+counter).mu_stat = 0;
                    h(start+counter).erp = EEG.data(channel, :, masked(jk,7));
                    h(start+counter).epochN = masked(jk,7);
                    h(start+counter).fullEp = squeeze(EEG.data(:, :, masked(jk,7)));
                end
            end
            save('dataframe.mat', 'h')
        end


        for i = 1:length(est_vec)

            % find in which epochs corresopond to condition A, which correspond to B:
            % [txt CWred CWnored SYL CLred CLnored]
            switch data_sep
                case 'median'
                    cutoff1 = nanmedian(masked(:, est_vec(i)));
                    cutoff2 = cutoff1;
                case 'mean'
                    cutoff1 = nanmean(masked(:, est_vec(i)));
                    cutoff2 = cutoff1;
                case 'quant'
                    med = nanmedian(masked(:, est_vec(i)));
                    cutoff1 = quantile(masked(:, est_vec(i)), 0.5);
                    cutoff2 = quantile(masked(:, est_vec(i)), 0.5);
            end
            
            % find epoch numbers above and below the cutoff value of WM load:
            epochNumbersInA = masked(find(masked(:, est_vec(i)) < cutoff1), 7); %EST{est_vec(i), 2});
            epochNumbersInB = masked(find(masked(:, est_vec(i))>= cutoff2), 7); %EST{est_vec(i), 2});

            % plot the ERPs:
            mFigure = subplot(3,2,i+Dspl);
            
            % get ERP data corresponding to epochs with high or low WM:
            Cond1 = squeeze(EEG.data(channel,:,epochNumbersInA));
            Cond2 = squeeze(EEG.data(channel,:,epochNumbersInB));
            save([filepath '/' 'ERP.mat'], 'Cond1','Cond2', 'epochNumbersInA', 'epochNumbersInB')
            
            % IMPLEMENT GA-BASED CLEANING:
            if genetic == true
                [epochNumbersInA, epochNumbersInB] = GA_SUPER(filepath);
                Cond1prime = squeeze(EEG.data(channel,:,epochNumbersInA));
                Cond2prime = squeeze(EEG.data(channel,:,epochNumbersInB));

                C1 = mean(Cond1,2); plot(EEG.times, C1); hold on; grid on
                C1 = mean(Cond1prime,2); plot(EEG.times, C1)         
                C2 = mean(Cond2,2); plot(EEG.times, C2)
                C2 = mean(Cond2prime,2); plot(EEG.times, C2)
            else
                C1 = mean(Cond1,2); plot(EEG.times, C1); hold on; grid on
                C2 = mean(Cond2,2); plot(EEG.times, C2)
            end
            % END of GA-BASED CLEANING
            
            % if AV.mat is not loaded yet, load:
            if exist('AV')==0
                load('/Users/RomanKoshkin/Documents/MATLAB/EEG/AV.mat')
            end
            
            % save data for later plotting grand average ERPs using the
            % legacy grand_avg_SUPER.m function.
            AV.(subject).(language).(EST{est_vec(i),1}).C1.data = squeeze(EEG.data(channel,:,epochNumbersInA));
            AV.(subject).(language).(EST{est_vec(i),1}).C1.labels = X1(epochNumbersInA,est_vec(i));
            AV.(subject).(language).(EST{est_vec(i),1}).C2.data = squeeze(EEG.data(channel,:,epochNumbersInB));      
            AV.(subject).(language).(EST{est_vec(i),1}).C2.labels = X1(epochNumbersInB,est_vec(i));
            
            % plot the ERPs for the current subject:
            titstr = [subject,' ', EST{est_vec(i),1},', ' num2str(data_sep)  ' cutoff: ' num2str(cutoff1) '-' num2str(cutoff2) ' / ', language];
                                                                    %EST{est_vec(i), 2}
            title(titstr, 'FontSize', 14);
            legend( strcat('low WM (', num2str(length(epochNumbersInA)),' trials)'),...
                    strcat('high WM (', num2str(length(epochNumbersInB)), ' trials)'));
            Avg = squeeze(EEG.data(channel,:,[epochNumbersInA; epochNumbersInB]));
            plot(EEG.times, mean(Avg,2), 'k', 'LineWidth', 2)
            line([0 0],[-2 2], 'LineWidth', 1.5);
            line([EEG.times(1) EEG.times(end)],[0 0], 'LineWidth', 1.5);

            % compute boundaries:
            Lb = round(-EEG.xmin*EEG.srate+LB*EEG.srate);
            Ub = round(-EEG.xmin*EEG.srate+UB*EEG.srate);
            
            % calculate condition-specific mean voltages within the window of interest:
            Cond1 = mean(Cond1(Lb:Ub,:),1);
            Cond2 = mean(Cond2(Lb:Ub,:),1);
            M1 = mean(Cond1);
            M2 = mean(Cond2);
            
            % run a independent samples t-test (exploratory):
            [h,p] = ttest2(Cond1,Cond2, 'Vartype', 'unequal');
            
            % add legends, make the plot pretty:
            str = strcat('MeanLow=', num2str(M1),', MeanHigh=', num2str(M2),char(10), ', Pval=', num2str(p));
            ylim=get(gca,'ylim');
            xlim=get(gca,'xlim');
            text(xlim(1)+5,ylim(2)-0.2,str)
            
            % highlight the window of interest:
            S.Vertices = [LB*1000 -2; LB*1000 2; UB*1000 2; UB*1000 -2];
            S.Faces = [1 2 3 4];
            S.EdgeColor = 'none';
            S.FaceAlpha = 0.25;
            patch(S)
            
            % clean the workspace:
            clear mFigure mTextBox temp
        end
    end
    
    % save
    save('/Users/RomanKoshkin/Documents/MATLAB/EEG/AV.mat', 'AV')
end