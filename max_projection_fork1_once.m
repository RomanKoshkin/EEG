%% see description in the e-mail, or at the bottom of the script
tic
clearvars -except h VV1
%% parameters:
subjects = {'KOK','KOS', 'ROM', 'POG', 'SHE', 'BUL', 'KOZ', 'ELT'} %, 'GRU'} % 
estimator = 'CLred'
langOfInterest = {'er', 'RE'}
language = {'er', 'RE'};
lwr = 0 % cut off the lower tail of the WM load distributions (set 0 if you don't want to)

q1 = 0.10+lwr; % lower quantile 50 for early topographies, 33 for P3b
q2 = 0.9; % upper quantile

%LB = [50  100 150 300];
%UB = [100 150 200 350];


LB = [0 100 180 320]; % in miliseconds
UB = [76 176 256 396];

% select the method implementation:
test = 1; % 1 - as in the paper, 0 - as you described in the email.
% regularization parameter:
lambda = 0.05;
% number of resamples for the randomization test:
perm = 2; 
%% internals:

% set counter to zero
uuu = 1;
P_rand = [0 0 -100 0 0];

% epoch times
t = [-100:4:496];

% load channel locations:
load ('chanlocs.mat');
channel = 22; %find(ismember({EEG.chanlocs.labels},'Cz') == 1);

if ~exist('h')
    load('dataframe_0.25-30Hz_500ms.mat') % load('dataframe_interp1_45Hz.mat')
    h(1:2) = [];
end
lenERP = length(h(2).erp);

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
        idx_low = find(idx & WW <= cutoff(1) & WW >= quantile(WW(idx), lwr));
        idx_high = find(idx & WW > cutoff(2));
        idx_mid = find(idx & WW <= cutoff(2) & WW > cutoff(1));
        [h(idx_low).tmplab] = deal({['low']});  % subjects{i} language{j}]});
        [h(idx_high).tmplab]= deal({['high']}); % subjects{i} language{j}]});
        [h(idx_mid).tmplab]= deal({['mid']});
    end
    end

% find epochs that match the conditions (low & high, E->R interpretation):
  % find epochs that match the conditions (low & high, E->R interpretation):
IDlow = find(ismember([h.tmplab], 'low') & ismember({h.lang}, langOfInterest));% & ismember({h.subj}, subjects));
IDhigh = find(ismember([h.tmplab], 'high') & ismember({h.lang}, langOfInterest));% & ismember({h.subj}, subjects));
% average these epochs:
low = mean(cat(3, h(IDlow).fullEp), 3);
high = mean(cat(3, h(IDhigh).fullEp), 3);

%% cleaning (remove epochs with voltage greater than +/- 50 uV)
% figure
% plot (mean(squeeze(low(22,:,:)),2)); hold on; plot(mean(squeeze(high(22,:,:)),2))
% 
% low = cat(3, h(IDlow).fullEp);
% xx = squeeze((low(22,:,:)));
% badLow = logical(sum(abs(xx) > 50));
% xx = xx(:, ~badLow);
% 
% high = cat(3, h(IDhigh).fullEp);
% yy = squeeze((high(22,:,:)));
% badHigh = logical(sum(abs(yy) > 30));
% yy = yy(:, ~badHigh);
% 
% figure
% plot(mean(xx,2)); hold on; plot(mean(yy,2))
%%



% computer the difference wave: 
DiffERP = low-high;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:4
% set the TR lower and upper bound, in samples:
[~, Lb] = min(abs(t-LB(k)));
[~, Ub] = min(abs(t-UB(k)));

% define TR and FR:
TR = [Lb:Ub];
FR = find(~ismember(1:lenERP, TR));

% correlation matrices for the target and flanker ranges:
C1 = DiffERP(:,TR)*DiffERP(:,TR)';
C2 = DiffERP(:,FR)*DiffERP(:,FR)';

% plot the ERPs at the channels of interest:
figure
subplot(1,3,1); q1 = gca; 
hold on; grid on
plot(q1, t, high(channel,:), 'r', 'LineWidth', 0.5)
plot(q1, t, low(channel,:), 'b', 'LineWidth', 0.5)
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

if test == 1 % DO IT THE WAY DESCRIBED IN THE PAPER (GEP)
            % Tikhonov regularized versions of the correlation matrices:
            C1 = C1 + lambda*trace(C1)/size(C1,1)*eye(size(C1,1));
            C2 = C2 + lambda*trace(C2)/size(C2,1)*eye(size(C2,1));
            % to find how to best project this data to maximize its variance (power),
            % we find the biggest eigenvector of its covariance matrix and project teh
            % data onto this eigenvector
            % v1 = v1*diag(1./sqrt(sum(v1.^2,1))); % This was in Sasha's script. What does this line do?
            [v1, d1] = eig(C1, C2);
            H = eye(length(v1)); % for this implementation we do not need whitening, so set it to I.

else        % DO IT ALEX'S WAY (AS IN THE EMAIL)
            H = inv(sqrtm(C2)); % get a whitening matrix based on the flanker data range (control)
            % H = H + lambda*trace(H)/size(H,1)*eye(size(H,1));
            y1 = H*DiffERP(:,TR); % put our TR data into this new space (basis)
            y2 = H*DiffERP(:,FR);
            Ry1 = y1*y1';   % get the covariance of this new data
            
            % Tikhonov regularization:
            Ry1 = Ry1 + lambda*trace(Ry1)/size(Ry1,1)*eye(size(Ry1,1));
            [v1, d1] = eig(Ry1);
end

[~, order] = sort(diag(d1),'descend');
v1 = v1(:,order);
d1 = d1(:,order);
% v1 = v1*diag(1./sqrt(sum(v1.^2,1)));

% plot the topographies:
subplot(1,3,3)
q3 = gca;
V = inv(v1); % topographies

xx = sign(mean(V(1,:)));
topoplot(V(1,:).*xx,chanlocs,'style','both','electrodes','labelpoint');
title (['Activation topographies, window: ' num2str([min(TR*4-100) max(TR*4-100)])])
q3.FontSize = 14;

subplot (1,3,2); q2 = gca;
plot(q2, t, v1(:,1)' .*xx * H*DiffERP, 'k', 'LineWidth', 2) %plot H-transformed ans projected data (difference)
hold on; grid on
plot(q2, t , v1(:,1)' .* xx * H*high, 'r', 'LineWidth', 0.5) %plot H-transformed ans projected data (highWM load)
plot(q2, t , v1(:,1)' .* xx * H*low, 'b', 'LineWidth', 0.5) %plot H-transformed ans projected data (lowWM load)
title('After projection')
q2.FontSize = 14; q2.YLim = [-3 3];

% highlight the TR

title('ERPs projected with computed spatial filters')

drawnow

% error('break')
%% permutation test

All = [IDlow, IDhigh]; % create a pile of all trial numbers in the conditions;

ss_true = v1(:,1)' * H * DiffERP;
true_score = abs(mean(ss_true(TR))-mean(ss_true(FR)));
lenA = length(IDhigh);
lenB = length(IDlow);
perm_score = zeros(perm, 1);
for j = 1:perm
    
    % take a random sample out of the pile:
    IDhigh_perm = randsample(All,lenA, 'false');
    IDlow_perm = All(~ismember(All, IDhigh_perm)); % randsample(All,lenB, 'false');
    
    % average trials within the appropriate conditions:
    low_perm = mean(cat(3, h(IDlow_perm).fullEp), 3);
    high_perm = mean(cat(3, h(IDhigh_perm).fullEp), 3);
    
    % get the difference wave for the resampled conditions: 
    DiffERP_perm = low_perm-high_perm;
    
    % correlation matrices for the target and flanker ranges:
    C1_perm = DiffERP_perm(:,TR)*DiffERP_perm(:,TR)';
    C2_perm = DiffERP_perm(:,FR)*DiffERP_perm(:,FR)';
    
    if test == 1 % DO IT THE WAY DESCRIBED IN THE PAPER (GEP)
            % Tikhonov regularized versions of the correlation matrices:
            C1_perm = C1_perm + lambda*trace(C1_perm)/size(C1_perm,1)*eye(size(C1_perm,1));
            C2_perm = C2_perm + lambda*trace(C2_perm)/size(C2_perm,1)*eye(size(C2_perm,1));
            % to find how to best project this data to maximize its variance (power),
            % we find the biggest eigenvector of its covariance matrix and project teh
            % data onto this eigenvector
            % v1 = v1*diag(1./sqrt(sum(v1.^2,1))); % This was in Sasha's script. What does this line do?
            [v1, d1] = eig(C1_perm, C2_perm);
            H = eye(length(v1)); % for this implementation we do not need whitening, so set it to I.

    else    % DO IT ALEX'S WAY (AS IN THE EMAIL)
            H = inv(sqrtm(C2)); % get a whitening matrix based on the flanker data range (control)
            y1 = H*DiffERP_perm(:,TR); % put our TR data into this new space (basis)
            y2 = H*DiffERP_perm(:,FR);
            Ry1 = y1*y1';   % get the covariance of this new data
            
            % Tikhonov regularization:
            Ry1 = Ry1 + lambda*trace(Ry1)/size(Ry1,1)*eye(size(Ry1,1));
            [v1, d1] = eig(Ry1);
    end

    % sort the eigenvectors and eigenvalues:
    [~, order] = sort(diag(d1),'descend');
    v1 = v1(:,order);
    d1 = d1(:,order);
    
    % project the permuted conditions onto our v1:
    ss_perm = v1(:,1)' * H * DiffERP_perm; % !!! i have removed H before DiffERP_perm
    
    % compute the difference change in the TR relative to FR:
    % perm_score (j) = abs(mean(ss_perm(TR))-mean(ss_perm(FR)));
    perm_score (j) = abs(mean(ss_perm(TR))) - abs(mean(ss_perm(FR)));
    
    % report progress:
    [num2str(j) ' of ' num2str(perm) ' samples complete']
end

% report a bootstrapped estimate of the p-value:
uuu = uuu+1; % increment the counter
P_rand(uuu,1:5) = [t(Lb) t(Ub) mean([t(Lb) t(Ub)])  mean([Lb Ub]) sum(perm_score > true_score)/perm]

toc
clf
VV(k,:,:) = V;
end

subplot(1,4,1)
topoplot(VV(1, 1,:).*-1,chanlocs,'style','both','electrodes','labelpoint');
title(['P1 (' num2str([LB(1) UB(1)]) ' ms)'], 'FontSize', 14)
subplot(1,4,2)
topoplot(VV(2, 1,:).*xx,chanlocs,'style','both','electrodes','labelpoint');
title(['N1 (' num2str([LB(2) UB(2)]) ' ms)'], 'FontSize', 14)
subplot(1,4,3)
topoplot(VV(3, 1,:).*xx,chanlocs,'style','both','electrodes','labelpoint');
title(['P2 (' num2str([LB(3) UB(3)]) ' ms)'], 'FontSize', 14)
subplot(1,4,4)
topoplot(VV(4, 1,:).*-1*xx,chanlocs,'style','both','electrodes','labelpoint');
title(['P3b (' num2str([LB(4) UB(4)]) ' ms)'], 'FontSize', 14)

%% e-mail:
% of course there is.
% the generalized eignevaue problem arises from the task to maximize Rayleigh quotient 
% Q = w'R_1w/w'R_2w aiming at finding the projection w which maximizes the power for 
% one data while minimizing it for the other. The data are reduced to their
% correlation matrices because of the power of w'X = w'XX'w = w'R_xw;
% Now, you are set to maximize this ratio and one way to solve this is to 
% work in the space where any direction w will give the same power 
% for the dataset that is in the denominator. 
% Then, in this space you simply need to find the direction that maximizes 
% the power for the first dataset(numerator)
% 
% to do that
% 1. Create a whiteninig transform(you can read about it leswhere) for the 
% second data H = inv(sqrtm(R_2)) which is inverse matrix square root 
% 2. take your data into this new space as 
% y1 = Hx1 and y2 = Hx2;
% 3. By contruction in this space the denominator of the Rayleigh quotient 
% will always be the same regardless of w (as long as the w'w = norm(w) = 1) 
% since w'Hxx'H'w = wHR_1H'w = w' I w = w'w = norm(w) = 1.  
% So, all you are  left with is to find max direction for 
% Ry1  = y1*y1' = HR_1H'  = inv(sqrtm(R_2))R_1inv(sqrtm(R_2)) 
% which is done via eigenvaue decomp of this matrix.  
% 
% QED