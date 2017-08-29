subjects = {'SHE', 'POG', 'BUL', 'KOS', 'ROM'};
for i = 1:length(subjects)
subject = subjects{i}
% 1 - CW with RED
% 4 - SYL
% 5 - CW without RED
% 6 - CL without RED
% 7 - CL with RED
% 8 - CW with RED (weighted)
% 9 - SYL (weighted)
% 10 - CW without RED (weighted)

estimate = 8; % 1, 4, 5 or 6 for unweighted CW, unweighted SYL, weighted CW and weighted SYL
path = ['/Users/RomanKoshkin/Documents/MATLAB/EEG/', subject, '_EEG/'];
load([path 'TimeCodeEVSdata_' subject '.mat'])
m = figure;
m.Name = subject;
s.er = [];
s.RE = [];
s.Mean.er = [];
s.SEM.er = [];
s.Mean.RE = [];
s.SEM.RE = [];
if estimate == 1 | estimate == 8
    hLim = 30;
else
    hLim = 60;
end
for i = 1:8
    a = s.(['text',num2str(i)]);
    b(:,1) = a(:,3) - a(:,2); % get chunk length
    a(:,8) = a(:,1).*b(:,1);  % compute weighed CW load in the chunk (with red)
    a(:,9) = a(:,4).*b(:,1);  % compute weighed SYL load in the chunk
    a(:,10) = a(:,5).*b(:,1);  % compute weighed CW load in the chunk (without red)  
    if strcmp(s.TextLabels(i,2),'er')
        s.er = [s.er; a(:,estimate)];      
    else
        s.RE = [s.RE; a(:,estimate)];
    end
    
    s.(['text',num2str(i)]) = a;
    clear a b
    subplot(6,2,i)
    x = s.(['text',num2str(i)])(:,estimate);
    y = hist(x,0:40)*100/length(x);
    bar(y)
    set(gca,'xlim',[0 hLim])
    set(gca,'ylim',[0 30])
    hold on
    plot ([mean(x) mean(x)], [0 30],'k-', 'LineWidth', 1)
    
    if strcmp(s.TextLabels(i,2),'er')
        s.Mean.er(end+1) = mean(x);
        s.SEM.er(end+1) = std(x)/sqrt(length(x));
    else
        s.Mean.RE(end+1) = mean(x);
        s.SEM.RE(end+1) = std(x)/sqrt(length(x));
    end
    
    ylim=get(gca,'ylim')-3;
    xlim=get(gca,'xlim');
    text(xlim(1),ylim(2),num2str(mean(x)))
    title([subject, ' ', num2str(i), ' ', char(s.TextLabels(i,1)), ' ', char(s.TextLabels(i,2))])
    grid on
end

for i = 1:2
    subplot(6,2,8+i)
    a = s.TextLabels{i,2};
    y = hist(s.(a), 0:40)*100/length(s.(a));
    bar(y)
    set(gca,'ylim',[0 30])
    set(gca,'xlim',[0 hLim])
    ylim=get(gca,'ylim')-3;
    xlim=get(gca,'xlim');
    x = mean(s.(a));
    text(xlim(1),ylim(2),num2str(x))
    hold on
    plot ([x x], [0 30],'k-', 'LineWidth', 1)
    title(['average ' a])
    grid on
end

for i = 1:2
    subplot(6,2,10 + i)
    a = s.TextLabels{i,2};
    errorbar([1:4],s.Mean.(a),s.SEM.(a),'-s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor','red')
    hold on
    plot ([1 4], [mean(s.Mean.(a)) mean(s.Mean.(a))],'--', 'LineWidth', 2)
    ax = gca;
    ax.XTick = ([1:4]);
    ax.XTickLabel = find(ismember(s.TextLabels(:,2), a));
    set(gca,'ylim',[1 4])
    title(['WM load dynamics by text ' a])
    xlabel('text number')
    grid on
end

clear a hLim m ax x y i xlim ylim estimate 
clc
end