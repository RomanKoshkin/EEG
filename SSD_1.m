
clear X z f Xc P Pc
ALLEEG = pop_delset(ALLEEG, 2:length(ALLEEG));
eeglab redraw
Fs = 250;
X = double(EEG.data (1:30,:));                  % get data
 
[de_srt, v_srt, topo_srt] = ssd(X,Fs,[8 12]); 
% search spatial filter to emphasize 8-12 band component 
% de_srt - ratio of component power to that of the entire data
% v_srt matrix of the spatial filter so that de_srt (1) = norm(v_srt(:,1)' * X)^2
% topo_srt - inverse of v_srt)

z = v_srt(:,1)'*X;              % applied the spatial filter (zero out the biggest 'component'
f = linspace(0,Fs,size(X,2));  % vector of freqs

load ('../EEG/chanlocs.mat')
chan = 'Pz';
channel = find(ismember({chanlocs.labels}, chan) == 1);

% figure
% plot(f,abs(fft(z)))                               % plot power spectrum
% plot(f,abs(fft(X(channel,:))), 'b')             % plot power spectrum
% Xc = X-topo_srt(:,1)*v_srt(:,1)'*X; % removing the biggest alpha-bearing comp. 

Xc = X-...
     topo_srt(:,1)*v_srt(:,1)'*X -...
     topo_srt(:,2)*v_srt(:,2)'*X -...; % removing 
     topo_srt(:,3)*v_srt(:,3)'*X;
 
% figure                                            %two most alpha-bearing components
% plot(f,abs(fft(X(1,:))),  'b')
% hold on 
% plot(f,abs(fft(Xc(1,:))), 'r')                             % plot power spectrum
% hold off

P = sum(X.*X,2);        % power over each electrode (dirty)
Pc = sum(Xc.*Xc,2);     % power over each electrode (clean)

figure; 
subplot(1,3,1)
topoplot(Pc,EEG.chanlocs,'style','both','electrodes','labelpoint'); title ('clean')
subplot(1,3,2)
topoplot(P,EEG.chanlocs,'style','both','electrodes','labelpoint'); title ('dirty')
subplot(1,3,3)
plot([1:30],P,'b', [1:30],Pc,'r')
legend ('dirty','clean')

figure
subplot(1,2,1)
spectopo(X(:,:), 0, EEG.srate,...
    'percent', 5,...
    'freqrange', [2 20],...
    'freq', [2,4,6,8,10,12],...
    'verbose', 'off',...
    'winsize', 150,...
    'chanlocs', EEG.chanlocs,...
    'nfft', 512,...
    'boundaries', [0,1, length(EEG.data)]);
grid on
title('DIRTY')

subplot(1,2,2)
spectopo(Xc(:,:), 0, EEG.srate,...
    'percent', 5,...
    'freqrange', [2 20],...
    'freq', [2,4,6,8,10,12],...
    'verbose', 'off',...
    'winsize', 150,...
    'chanlocs', EEG.chanlocs,...
    'nfft', 512,...
    'boundaries', [0,1, length(EEG.data)]);
grid on
title('CLEAN')

% save data:

EEG.data = Xc;
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', 'Clean');
EPOCH_SUPER
cl = mean(EEG.data(22,:,:),3);

EEG = eeg_retrieve(ALLEEG, 1); CURRENTSET = 1; eeglab redraw
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', 'Dirty');
stri = strfind({ALLEEG.setname}, 'Dirty');
addre = find(~cellfun(@isempty,stri));
EEG = eeg_retrieve(ALLEEG, addre); CURRENTSET = addre; eeglab redraw

EPOCH_SUPER
dirty = mean(EEG.data(22,:,:),3);

plot(dirty);hold on; plot(cl)