%% Run this after NEW_GRAND_AVERAGE.m

x = zeros(30,150);
for i = a;                   % LOW
    x = x + h(i).fullEp;
end
x = x/length(a);
% figure; plot(mean(x,1))
fig = figure;

% j = [-0.5 -0.3 -0.1 -0.1 -0.1 -0.1 -0.1];
t = [-100:4:496]';
time = [0 25 50 75 100 125 150 175 200;]
times = [26 32 39 45 51 57 64 70 76];
for i = 1:length(times)
    subplot (2,length(times),i)
    topoplot(x(:,times(i)), EEG.chanlocs, 'maplimits',[-0.8 0.8]);
    title([num2str(time(i)) ' ms'], 'FontSize', 14)
end
%%
x = zeros(30,150);
for i = b;                   % High
    x = x + h(i).fullEp;
end
x = x/length(b);

t = [-100:4:496]';
times = [26 32 39 45 51 57 64 70 76];
for i = 1:length(times)
    subplot (2,length(times),i + length(times))
    topoplot(x(:,times(i)), EEG.chanlocs, 'maplimits',[-0.8 0.8]);
    title([num2str(time(i)) ' ms'], 'FontSize', 14)
end

cbar('vert',0,[-1 1])
fig.Position = [0 0 1440 284];
fig.PaperSize = [17 3.9];
fig.PaperType = 'usletter';
fig.PaperUnits = 'inches';
% print (fig, '-r600', '-dpdf',  '-r0')