% this script generates (through the function generate_new.m) synthetic EEG
% with 59 channels. The 32-channel electrodes are 1:28. The sampling rate 500 Hz.
% excluding FT9 and FT10 (they were references)
clear
%% parameters:
testprojected = 1; % 0 if you want to do a randomization test on non-projected ERPs
channel = 21;
component_of_interest = 2; % will projecting on the eigenvector, corresponding the the biggest eigenvalue
win = 40; % TR window in samples
step = 4; % scanning the entire ERP waveform in steps of 4 samples
% select the method implementation:
test = 1; % 1 - as in the paper, 0 - as you described in the email.
% regularization parameter:
lambda = 0.05; 
% number of resamples for the randomization test:
perm = 75; 
%% internals:

% generate synthetic data:
epochs = 100; % how many epochs to generate
LAMBDA = 100; % noise parameter 
lenERP = 500; % length of the ERP waveform (cannot set)
[Data_s, Data_d, Ns, T] = generate_new(0, epochs, LAMBDA, 3);
Data_s = Data_s(1:28,:);
Data_d = Data_d(1:28,:);
load('chanlocs_28.mat')

j = 0;
dev = zeros(28,500,epochs);
sta = zeros(28,500,epochs);
for i = 1:500:(length(Data_s-499))
    j = j + 1;
    dev(:,:,j) = Data_d(:,i:(i+499));
    sta(:,:,j) = Data_s(:,i:(i+499));
end

% set counter to zero
uuu = 1;
P_rand = [0 0 -100 0 0];

% epoch times
t = [1:500];

% load channel locations:
load ('chanlocs_28.mat');

low = mean(dev,3);
high = mean(sta,3);

% computer the difference wave: 
DiffERP = low-high;

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
for u = 0:step:(lenERP-win)      
    Lb = 1 + u;
    Ub = win + u;
    clf

    % define TR and FR:
    TR = [Lb:Ub];
    FR = find(~ismember(1:lenERP, TR));

    % correlation matrices for the target and flanker ranges:
    C1 = DiffERP(:,TR)*DiffERP(:,TR)';
    C2 = DiffERP(:,FR)*DiffERP(:,FR)';

    % plot the ERPs at the channels of interest:
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
    q1.FontSize = 14; q1.YLim = [-1.1 1.1];

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
topoplot(V(component_of_interest,:).*xx,chanlocs,'style','both','electrodes','labelpoint');
title 'Activation topographies'
q3.FontSize = 14;

subplot (1,3,2); q2 = gca;
plot(q2, t, v1(:,component_of_interest)' .*xx * H*DiffERP, 'k', 'LineWidth', 2) %plot H-transformed ans projected data (difference)
hold on; grid on
plot(q2, t , v1(:,component_of_interest)' .* xx * H*high, 'r', 'LineWidth', 0.5) %plot H-transformed ans projected data (highWM load)
plot(q2, t , v1(:,component_of_interest)' .* xx * H*low, 'b', 'LineWidth', 0.5) %plot H-transformed ans projected data (lowWM load)
title('After projection')
q2.FontSize = 14; q2.YLim = [-1.1 1.1];

% highlight the TR

title('ERPs projected with computed spatial filters')

drawnow
% error('break')
%% permutation test
    allData = cat(3, dev, sta);
    IDdev = [1:epochs];
    IDsta = [epochs+1:epochs*2];
    All = [IDdev, IDsta]; % create a pile of all trial numbers in the conditions;

    ss_true = v1(:,component_of_interest)' * H * DiffERP;
    true_score = abs(mean(ss_true(TR))-mean(ss_true(FR)));
    lendev = length(IDdev);
    lensta = length(IDsta);
    perm_score = zeros(perm, 1);
    for j = 1:perm

        % take a random sample out of the pile:
        IDdev_perm = randsample(All,lendev, 'false');
        IDsta_perm  = randsample(All,lensta, 'false');

        % average trials within the appropriate conditions:
        dev_perm = squeeze(mean(allData(:,:,IDdev_perm),3));
        sta_perm = squeeze(mean(allData(:,:,IDsta_perm),3));

        % get the difference wave for the resampled conditions: 
        DiffERP_perm = dev_perm-sta_perm;

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
        
        if testprojected == 0
            v1 = ones(size(v1));
        end

        % project the permuted conditions onto our v1:
        ss_perm = v1(:,component_of_interest)' * H * DiffERP_perm; % !!! i have removed H before DiffERP_perm

        % compute the difference change in the TR relative to FR:
        % perm_score (j) = abs(mean(ss_perm(TR))-mean(ss_perm(FR)));
        perm_score (j) = abs(mean(ss_perm(TR))) - abs(mean(ss_perm(FR)));

        % report progress:
        % [num2str(j) ' of ' num2str(perm) ' samples complete']
    end

    % report a bootstrapped estimate of the p-value:
    uuu = uuu+1; % increment the counter
    P_rand(uuu,1:5) = [t(Lb) t(Ub) mean([t(Lb) t(Ub)])  mean([Lb Ub]) sum(perm_score > true_score)/perm];
    toc
end
%% plot topographies:

% plot the stem plot for p-values:
figure
stem(P_rand(:,3),P_rand(:,5))
hold on
plot([-100 500], [0.01,0.01])
hold on
stem(P_rand(find(P_rand(:,5)<=0.01),3), P_rand(find(P_rand(:,5)<=0.01),5))