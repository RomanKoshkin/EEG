% This script performs a randomization test of equal latencies in the
% windows of interest. The null hypothesis is that the latencies are the
% same

clear
tic
dataframe = 'dataframe_0.25-30Hz_500ms.mat';
load(dataframe)
toc
%% parameters:
subjects = {'KOK', 'KOS', 'ROM', 'POG', 'ELT', 'SHE', 'BUL', 'KOZ', 'GRU'}
t = [-100:4:496];
estimators = {'CLred', 'CWred','CLnored','CWnored', 'SYL'};
window = 54:73; % 54:73 (doesn't give bad negative diff in the N1)
% ;%37:44; %58:64 for N1 %57:67 % for P1 REPLACE ALL min WITH MIN, AND VICE VERSA

% number of permutations:
perm = 1000;

counter = 0;
p_val = struct;

%%
for kkk = 1:length(estimators)
    clearvars -except subjects estimators kkk counter p_val h window t perm
    estimator = estimators{kkk};
    language = {'er', 'RE'};
    languages = {'er', 'RE', 'both'};
    q1 = 0.1; % lower quantile
    q2 = 0.9; % upper quantile

    %% internals:

    lenERP = length(h(2).erp);
    
    % set counter to zero
    uuu = 1;
    
    % initialize the first row of the matrix of p-values:
    P_rand = [0 0 -100 0 0];

    % load channel locations:
    load ('chanlocs.mat');
    channel = 22; %find(ismember({EEG.chanlocs.labels},'Cz') == 1);

    % NOW COMPUTE THE SUBJECT- AND LANGUAGE-SPECIFIC QUANTILES AND CUTOFFS:    
    y = {};
    iter = 0;
    WW = [h.(estimator)];
    SSubj = {h.subj};
    Slang = {h.lang};
    [h.tmplab] = deal({'none'});
    for i = 1:length(subjects)
        id2 = ismember(SSubj, subjects{i});
        for j=1:length(language)
            id3 = ismember(Slang, language{j});
            iter = iter + 1;
            idx = id2 & id3;

            y{iter, 1} = [subjects{i} language{j}];
            y{iter,2} = WW(idx); % mean voltage in ERP window subjXlang
            cutoff = quantile(y{iter, 2}, [q1 q2]);
            y{iter,3} = cutoff;
            idx_low = find(idx & WW <= cutoff(1));
            idx_high = find(idx & WW > cutoff(2));
            idx_mid = find(idx & WW <= cutoff(2) & WW > cutoff(1));
            [h(idx_low).tmplab] = deal({['low']});  % subjects{i} language{j}]});
            [h(idx_high).tmplab]= deal({['high']}); % subjects{i} language{j}]});
            [h(idx_mid).tmplab]= deal({['mid']});
        end
    end

    % find epochs IDs that match the conditions (low & high, E->R interpretation):
    allIDlow = find(ismember([h.tmplab], 'low'));
    allIDhigh = find(ismember([h.tmplab], 'high'));

    REIDlow = find(ismember([h.tmplab], 'low') & ismember({h.lang}, 'RE'));
    REIDhigh = find(ismember([h.tmplab], 'high') & ismember({h.lang}, 'RE'));

    erIDlow = find(ismember([h.tmplab], 'low') & ismember({h.lang}, 'er'));
    erIDhigh = find(ismember([h.tmplab], 'high') & ismember({h.lang}, 'er'));

    % average these epochs:
    alllow = mean(cat(1, h(allIDlow).erp), 1);
    allhigh = mean(cat(1, h(allIDhigh).erp), 1);

    erlow = mean(cat(1, h(erIDlow).erp), 1);
    erhigh = mean(cat(1, h(erIDhigh).erp), 1);
    figure; plot(t, erlow, 'b', t, erhigh, 'r')

    RElow = mean(cat(1, h(REIDlow).erp), 1);
    REhigh = mean(cat(1, h(REIDhigh).erp), 1);

    % find N1 peak latencies in the averaged epochs:
    truLo_all = t(find(alllow==min(alllow(window))));
    truHi_all = t(find(allhigh==min(allhigh(window))));

    truLo_er = t(find(erlow==min(erlow(window))));
    truHi_er = t(find(erhigh==min(erhigh(window))));

    truLo_RE = t(find(RElow==min(RElow(window))));
    truHi_RE = t(find(REhigh==min(REhigh(window))));
    peaks = [truLo_RE truHi_RE; truLo_er truHi_er; truLo_all truHi_all];

    ['tru_er: '     num2str(truLo_er)   '  ' num2str(truHi_er)]
    ['tru_RE: '     num2str(truLo_RE)   '  ' num2str(truHi_RE)]
    ['tru_all: '    num2str(truLo_all)  '  ' num2str(truHi_all)]

    %% Clean by removing epochs in which the voltage exceeds the threshold
    % co = 25; % voltage threshold
    % badID = sum(RElow < -co,2) & sum(RElow > co,2);
    % allID = 1:length(RElow);
    % [num2str(sum(badID)/length(allID)*100) ' % removed']
    % goodID = allID(~badID);
    % plot(mean(RElow(allID,:),1))
    % hold on
    % plot(mean(RElow(goodID,:),1))
    %% permutation test

    for LANG = 1:3

        switch(LANG)
            case 1
            IDlow = REIDlow;
            IDhigh = REIDhigh;
            true_score = truHi_RE-truLo_RE;
            case 2
            IDlow = erIDlow;
            IDhigh = erIDhigh;
            true_score = truHi_er-truLo_er;
            case 3
            IDlow = allIDlow;
            IDhigh = allIDhigh;
            true_score = truHi_all-truLo_all;
        end

        All = [IDlow, IDhigh]; % create a pile of all trial numbers in the conditions;

        lenA = length(IDhigh);
        lenB = length(IDlow);
        P_rand = zeros(1, perm);
        clear perm_score ID_perm IDlow_perm low_perm high_perm permLo permHi
        counter = counter + 1;
        for j = 1:perm

            % take a random sample out of the pile:
            IDhigh_perm = randsample(All,lenA, 'false');
            % IDlow_perm = randsample(All,lenB, 'false');
            IDlow_perm  = All(~ismember(All, IDhigh_perm));

            % average trials within the appropriate conditions:
            low_perm = mean(cat(1, h(IDlow_perm).erp), 1);
            high_perm = mean(cat(1, h(IDhigh_perm).erp), 1);

            % get the difference wave for the resampled conditions: 
            permLo = t(find(low_perm==min(low_perm(window))));
            permHi = t(find(high_perm==min(high_perm(window))));

            % compute the difference change in the TR relative to FR:
            perm_score (j) = abs(permHi - permLo);

            P_rand(j) = sum(perm_score > true_score)/j;

            % create and maximize a figure if desn't exist yet:
            if ~exist('WWW')
                WWW = figure;
                set(WWW,'units','normalized','outerposition',[0 0 1 1])
            end

            % plot the permuted data:
            subplot(1,2,1)
            plot(t, low_perm, t, high_perm)
            title([languages(LANG) '  ' estimator])
            ax = gca;
            ax.FontSize = 14        

            % display the patch:
            S.Vertices = [t(window(1)) -2; t(window(1)) 2; t(window(end)) 2; t(window(end)) -2];
            S.Faces = [1 2 3 4];
            S.EdgeColor = 'none';
            S.FaceAlpha = 0.25;
            patch(S)

            grid on
            subplot(1,2,2)
            plot(1:j, P_rand(1:j))
            drawnow

            % report progress:
        end
        % p_val(LANG) = P_rand(j);
        % title(['lang:' num2str(LANG)...
        %     '; p = ' num2str(P_rand(j))...
        %     '; perm.: ' num2str(perm)...
        %     '; quant: ' num2str(q1) '-' num2str(q2) ' '...
        %     estimator])
        counter
        p_val(counter).estimator = estimator;
        p_val(counter).P_rand = P_rand(j);
        p_val(counter).lang = languages{LANG};
        p_val(counter).quant = [num2str(q1), '_' num2str(q2)];
        p_val(counter).true_score = true_score;
        p_val(counter).times = [num2str(t(58)) '_' num2str(t(64))];
        p_val(counter).peaks = peaks(LANG,:);
    end
end
disp('Done!')