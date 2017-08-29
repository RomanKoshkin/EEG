%% see description in you e-mail, or at the bottom of the script
clearvars -except s
%% parameters:
if ~exist('s')
    load('/Users/RomanKoshkin/Downloads/dataframe.mat')
end
test = 0; % 1 - as in the paper, 0 - as you described in the email.
place = 0;

LB = 0.01; % seconds. Lower bound for measuring the mean votage in the N1 interval
UB = 0.03;
lambda = 0.1;
iter = 50; % how many times to resample in the test
%% 
% internals:
load ('chanlocs.mat');
xmin = -0.100;
srate = 250;
channel = 22; %find(ismember({EEG.chanlocs.labels},'Cz') == 1);
Lb = round(-xmin*srate+LB*srate);
Ub = round(-xmin*srate+UB*srate);
%%
% times
t = [-100:4:396];
TR = [Lb:Ub];
FR = [1:Lb-1 Ub+1:length(t)];

% average matrices:
IDlow = find(ismember([s.tmplab], 'low') & ismember({s.lang}, 'er'));
IDhigh = find(ismember([s.tmplab], 'high') & ismember({s.lang}, 'er'));

low = mean(cat(3, s(IDlow).fullEp), 3);
high = mean(cat(3, s(IDhigh).fullEp), 3);

% difference waves: 
DiffERP = low-high;

% correlation matrices for the target and flanker ranges:
C1 = DiffERP(:,TR)*DiffERP(:,TR)';
C2 = DiffERP(:,FR)*DiffERP(:,FR)';

subplot(1,3,1); q1 = gca; 
plot(t, DiffERP(channel,:), 'k')
hold on; grid on
plot(t, high(channel,:), 'r')
plot(t, low(channel,:), 'b')
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
            H = C2^(1/2); % get a whitening matrix based on the flanker data range (control)
            H = inv(H + lambda*trace(H)/size(H,1)*eye(size(H,1)));
            y1 = H*DiffERP(:,TR); % put our TR data into this new space (basis)
            y2 = H*DiffERP(:,FR);
            Ry1 = y1*y1';   % get the covariance of this new data
            Ry2 = y2*y2';
            % Tikhonov regularization:
            Ry1 = Ry1 + lambda*trace(Ry1)/size(Ry1,1)*eye(size(Ry1,1));
            Ry2 = Ry2 + lambda*trace(Ry2)/size(Ry2,1)*eye(size(Ry2,1));
            [v1, d1] = eig(Ry1);
end

[~, order] = sort(diag(d1),'descend');
v1 = v1(:,order);
d1 = d1(:,order);

subplot (1,3,2); q2 = gca;
plot(t ,v1(:,1)'* H*DiffERP, 'k') %plot H-transformed ans projected data (difference)
hold on; grid on
plot(t ,v1(:,1)'* H*high, 'r') %plot H-transformed ans projected data (highWM load)
plot(t ,v1(:,1)'* H*low, 'b') %plot H-transformed ans projected data (lowWM load)
title('After projection')
q2.FontSize = 14; q2.YLim = [-3 3];

S.Vertices = [t(Lb) -2; t(Lb) 2; t(Ub) 2; t(Ub) -2];
S.Faces = [1 2 3 4];
S.EdgeColor = 'none';
S.FaceAlpha = 0.25;
patch(S)
% set(gca,'YDir','reverse')
title('ERPs projected with computed spatial filters')

subplot(1,3,3)
q3 = gca;
V = inv(v1); % topographies
topoplot(V(1,:),chanlocs,'style','map','electrodes','labelpoint');
title 'Activation topographies'
q3.FontSize = 14;

drawnow
% error('break')
%% permutation test

All = [IDlow, IDhigh]; % all trial numbers in given conditions;

ss_true = v1(:,1)'* H * DiffERP;
true_score = abs(mean(ss_true(TR))-mean(ss_true(FR)));
lenA = length(IDhigh);
lenB = length(IDlow);
perm_score = zeros(iter, 1);
for j = 1:iter
    
    % take a random sample out of the pile:
    IDhigh = randsample(All,lenA, 'false');
    IDlow  = randsample(All,lenB, 'false');
    
    % average trials within the appropriate conditions:
    low = mean(cat(3, s(IDlow).fullEp), 3);
    high = mean(cat(3, s(IDhigh).fullEp), 3);
    
    % get the difference wave for the resampled conditions: 
    DiffERP = low-high;
    
    % correlation matrices for the target and flanker ranges:
    C1 = DiffERP(:,TR)*DiffERP(:,TR)';
    C2 = DiffERP(:,FR)*DiffERP(:,FR)';
    
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

    else    % DO IT ALEX'S WAY (AS IN THE EMAIL)
            H = C2^(1/2); % get a whitening matrix based on the flanker data range (control)
            H = inv(H + lambda*trace(H)/size(H,1)*eye(size(H,1)));
            y1 = H*DiffERP(:,TR); % put our TR data into this new space (basis)
            y2 = H*DiffERP(:,FR);
            Ry1 = y1*y1';   % get the covariance of this new data
            Ry2 = y2*y2';
            % Tikhonov regularization:
            Ry1 = Ry1 + lambda*trace(Ry1)/size(Ry1,1)*eye(size(Ry1,1));
            Ry2 = Ry2 + lambda*trace(Ry2)/size(Ry2,1)*eye(size(Ry2,1));
            [v1, d1] = eig(Ry1);
    end

    % sort the eigenvectors and eigenvalues:
    [~, order] = sort(diag(d1),'descend');
    v1 = v1(:,order);
    d1 = d1(:,order);
    
    
    % project the permuted conditions onto our v1:
    ss_perm = v1(:,1)'* H * DiffERP;
    
    % compute the difference change in the TR relative to FR:
    perm_score (j) = abs(mean(ss_perm(TR))-mean(ss_perm(FR)));
    
    % report progress:
    [num2str(j) ' of ' num2str(iter) ' samples complete']
end

% report a bootstrapped estimate of the p-value:
P_rand = sum(perm_score > true_score)/iter

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