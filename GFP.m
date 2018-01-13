%% global field power:
% https://sccn.ucsd.edu/pipermail/eeglablist/2008/002132.html
% 
figure;
plot(EEG.times, mean(EEG.data,3), 'k');
hold on;
plot(EEG.times, std(mean(EEG.data,3)), 'r', 'linewidth', 3);

%% explore channel spectra and topographies at set frequencies:
clc
spectopo(EEG.data(:,:), 0, EEG.srate,...
    'chanlocs', EEG.chanlocs,...
    'freq', [1 3 5 7 9 11 13 15 17 19 21],...
    'percent', 10,...
    'freqrange', [0 40],...
    'winsize', 150,...
    'nfft', 512,...
    'mapchans', 22);

%% compute GFP based on the h struct
for i = [1:6480, 9109:53029] % without GRU
    XX(:,:,i) = h(i).fullEp;
end

low_mrk = find(ismember([h.tmplab], 'low'));
high_mrk = find(ismember([h.tmplab], 'high'));
mid_mrk = find(ismember([h.tmplab], 'mid'));

low = mean(XX(:,:,low_mrk),3);
mid = mean(XX(:,:,mid_mrk),3);
high = mean(XX(:,:,high_mrk),3);

figure
plot([-100:4:496], low', 'k', 'linewidth', 1);
hold on
plot([-100:4:496], std(low), 'b', 'linewidth', 1);
plot([-100:4:496], std(mid), 'k', 'linewidth', 1);
plot([-100:4:496], std(high),'r', 'linewidth', 1);