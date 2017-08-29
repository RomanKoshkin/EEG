tic
clear
%% parameters:
subjects = {'KOK','KOS', 'ROM', 'POG', 'ELT', 'SHE', 'BUL', 'KOZ'} % 'GRU'
estimator = 'CLred'
langOfInterest = {'er', 'RE'}
language = {'er', 'RE'};

win = 10;
step = 1;

q1 = 0.10; % lower quantile
q2 = 0.90; % upper quantile

% number of resamples for the randomization test:
perm = 200; 
%% internals:

% load the dataset is not loaded yet:
if ~exist('h')
    load('dataframe_0.25-30Hz_500ms.mat')
    h(1:2) = [];
end
lenERP = length(h(2).erp);
% set counter to zero
uuu = 1;
P_rand = [0 0 -100 0 0];

% epoch times
t = [-100:4:496];

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
% find epochs that match the conditions (low & high, E->R interpretation):
IDlow = find(ismember([h.tmplab], 'low') & ismember({h.lang}, langOfInterest));% & ismember({h.subj}, subjects));
IDhigh = find(ismember([h.tmplab], 'high') & ismember({h.lang}, langOfInterest));% & ismember({h.subj}, subjects));

% average these epochs:
low = mean(cat(1, h(IDlow).erp), 1);
high = mean(cat(1, h(IDhigh).erp), 1);

% computer the difference wave: 
DiffERP = low-high;


tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% for u = 0:2:lenERP - win % we'll shift the TR in steps of 2 samples (8 ms)
for u = 0:step:(lenERP-win)      
    Lb = 1 + u;
    Ub = win + u;
    clf

    % define TR and FR:
    TR = [Lb:Ub];
    FR = find(~ismember(1:lenERP, TR));


    % plot the ERPs at the channels of interest:
    subplot(1,3,1); q1 = gca; 
    hold on; grid on
    plot(q1, t, high, 'r', 'LineWidth', 0.5)
    plot(q1, t, low, 'b', 'LineWidth', 0.5)
    if exist('P_rand')
        plot(q1, P_rand(:, 3), P_rand(:,5), 'k', 'LineWidth', 3)
    end
    S.Vertices = [t(Lb) -2; t(Lb) 2; t(Ub) 2; t(Ub) -2];
    S.Faces = [1 2 3 4];
    S.EdgeColor = 'none';
    S.FaceAlpha = 0.25;
    patch(S)
    title 'Before projection, Cz'
    q1.FontSize = 14; q1.YLim = [-2 2];


drawnow
% error('break')
    %% permutation test

    All = [IDlow, IDhigh]; % create a pile of all trial numbers in the conditions;

    ss_true = DiffERP;
    true_score = abs(mean(ss_true(TR))-mean(ss_true(FR)));
    lenA = length(IDhigh);
    lenB = length(IDlow);
    perm_score = zeros(perm, 1);
    for j = 1:perm

        % take a random sample out of the pile:
        IDhigh_perm = randsample(All,lenA, 'false');
        IDlow_perm = All(~ismember(All, IDhigh_perm)); % randsample(All,lenB, 'false');
        
        % average trials within the appropriate conditions:
        low_perm = mean(cat(1, h(IDlow_perm).erp), 1);
        high_perm = mean(cat(1, h(IDhigh_perm).erp), 1);

        % get the difference wave for the resampled conditions: 
        DiffERP_perm = low_perm-high_perm;

        % project the permuted conditions onto our v1:
        ss_perm = DiffERP_perm; % !!! i have removed H before DiffERP_perm

        % compute the difference change in the TR relative to FR:
        perm_score (j) = abs(mean(ss_perm(TR))) - abs(mean(ss_perm(FR)));

        % report progress:
        % [num2str(j) ' of ' num2str(perm) ' samples complete']
    end

    % report a bootstrapped estimate of the p-value:
    uuu = uuu+1; % increment the counter
    P_rand(uuu,1:5) = [t(Lb) t(Ub) mean([t(Lb) t(Ub)])  mean([Lb Ub]) sum(perm_score > true_score)/perm]
    toc
end

% plot the stem plot for p-values:

figure
stem(P_rand(:,3),P_rand(:,5))
hold on
plot([-100 500], [0.05,0.05])
hold on
stem(P_rand(find(P_rand(:,5)<=(0.05/perm)),3), P_rand(find(P_rand(:,5)<=(0.05/perm)),5))
title('RAW RANDOMIZATION P-VALUES')
% BENJAMINI-HOCHBERG CORRECTED P-VALUES:
% plot the stem plot for p-values:
P_rand_corr = readtable('/Users/RomanKoshkin/Documents/R/BH.csv');
P_rand_corr = sortrows(P_rand_corr, 'onset');
figure
stem(P_rand_corr.midpoint, P_rand_corr.BH)
hold on
plot([-100 500], [0.05,0.05])
hold on
stem(P_rand_corr.midpoint(find(P_rand_corr.BH<=0.05)),...
P_rand_corr.BH(find(P_rand_corr.BH<=0.05)))
title('BENJAMINI-HOCHBERG CORRECTED P-VALUES')