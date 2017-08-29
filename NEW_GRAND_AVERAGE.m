clearvars -except h
win = 'mu_N1';
kernel_len = 1;
% N1
LB_N1 = 0.13;   % 0.13
UB_N1 = 0.16;   % 0.15
% P1
LB_P1 = 0.05;   % 0.05
UB_P1 = 0.08;   % 0.07
Xtab = table;
Xtab.times = [-100:4:496]';
estimator = 'CLred';
subjects =              {'KOK', 'ELT', 'KOZ', 'POG', 'KOS', 'ROM', 'SHE', 'BUL','GRU'};
subjectsofinterest =    {'KOK', 'ELT', 'KOZ', 'POG', 'KOS', 'ROM', 'SHE', 'BUL','GRU'};
WMloadsofinterest = {'high', 'mid', 'low'};
language = {'er', 'RE'};
chan = 'Cz';
if ~exist('h')
    load('dataframe_0.25-30Hz_500ms.mat') % dataframe_20Hz_interp1.mat
    % dataframe_0.25-30Hz_500ms.mat - is the same as
    % dataframe_0.25-30Hz.mat, but the epoch length post-stim is 500 ms,
    % not 400 ms. Both are FIR filtered between 0.25 and 30 Hz.
    h(1:2) = [];
end
%% internals:
if isfield(h, 'tmplab')==0
    [h.tmplab] = deal({['a']});
end
load ('chanlocs.mat');
channel = find(ismember({chanlocs.labels}, chan) == 1);
prompt = {'Q1:','Q2:', 'step:', 'full sim?', 'plot ERPs?', 'updata df 4 R?:'};
dlg_title = 'Set parameters';
num_lines = 1;
defaultans = {'0.15', '0.85', '0.05', '0', '1', '0'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
if str2num(answer{4})==1
    qStep = str2num(answer{3});
    Q2 = str2num(answer{1}):qStep:str2num(answer{2});
    C = zeros(length(Q2),length(Q2), 2);
    E = zeros(length(Q2),length(Q2), 2);
    
else
    Q2 = str2num(answer{2});
    Q1 = str2num(answer{1});
end
saveR = str2num(answer{6});
BL = 0;
sw1 = 0;
bigLoop = length(Q2)^2/2;
plotting = str2num(answer{5});
for jj = 1:length(Q2)
    tic
    if exist('qStep')
        Q1 = [str2num(answer{1}):qStep:Q2(jj)];
    end
    for ii = 1:length(Q1)
        BL = BL+1;
        

q1 = Q1(ii); % lower quantile
q2 = Q2(jj);


Lb_N1 = round(0.1*250+LB_N1*250);
Ub_N1 = round(0.1*250+UB_N1*250);
Lb_P1 = round(0.1*250+LB_P1*250);
Ub_P1 = round(0.1*250+UB_P1*250); 
cutoffOVERALL = quantile([h.(estimator)], [q1 q2]);

% recompute the means within the window of interest:
    for i = 1:length(h)
        h(i).erp = h(i).fullEp(channel,:);
    end
    mu = reshape([h.erp], length(Xtab.times), []);
    tmp_N1 = num2cell(mean(mu(Lb_N1:Ub_N1,:),1));
    tmp_P1 = num2cell(mean(mu(Lb_P1:Ub_P1,:),1)); clear mu;
    [h.mu_N1] = deal(tmp_N1{:}); clear tmp_N1
    [h.mu_P1] = deal(tmp_P1{:}); clear tmp_P1
    
%     % check just in case the average computation:
%     plot(mid_mu); hold on;
%     MM = mean(mid_mu(Lb_N1:Ub_N1));
%     plot(1:150,linspace (MM, MM, 150))
%     plot(linspace (Lb_N1, Lb_N1, 150), linspace (-0.8, 0.8, 150))
%     plot(linspace (Ub_N1, Ub_N1, 150), linspace (-0.8, 0.8, 150))

% NOW COMPUTE THE SUBJECT- AND LANGUAGE-SPECIFIC QUANTILES AND CUTOFFS:  
%% alternative 1 (quantile based):
y = {};
iter = 0;
WW = [h.(estimator)];
SSubj = {h.subj};
Slang = {h.lang};
[h.tmplab] = deal({''});
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
%         histogram(WW(idx))
%         title([subjects{i} ' ' language{j}])
        idx_low = find(idx & WW <= cutoff(1));
        idx_high = find(idx & WW > cutoff(2));
        idx_mid = find(idx & WW <= cutoff(2) & WW > cutoff(1));
        [h(idx_low).tmplab] = deal({'low'});  % subjects{i} language{j}]});
        [h(idx_high).tmplab]= deal({'high'}); % subjects{i} language{j}]});
        [h(idx_mid).tmplab]= deal({'mid'});
    end
end
%% alternative 2 (average, quantile-free)
% % y = {};
% % iter = 0;
% % WW = [h.(estimator)];
% % SSubj = {h.subj};
% % Slang = {h.lang};
% % [h.tmplab] = deal({''});
% % 
% % cutoff = quantile(WW, [q1 q2]);
% % idx_low = find(WW < cutoff(1));
% % idx_high = find(WW >= cutoff(2));
% % idx_mid = find(WW >= cutoff(1) & WW < cutoff(2));
% % 
% % [h(idx_low).tmplab] = deal({'low'});
% % [h(idx_high).tmplab]= deal({'high'});
% % [h(idx_mid).tmplab]= deal({'mid'});

%%
g = 1;
    if plotting==1;
    if ~exist('WWW')
        WWW = figure;
        plot1 = subplot(2,3,1);
        plot2 = subplot(2,3,2);
        plot3 = subplot(2,3,3);
        plot4 = subplot(2,3,4);
        plot5 = subplot(2,3,5);
        plot6 = subplot(2,3,6);
    end
    if sw1 == 0
        boxplot (plot3, [h.(estimator)], {{h.subj}, {h.lang}});
        grid on
        txt11 = text(4.6483, 8.9443, 'WM load by subject and direction of translation');
        sw1 = 1;
    end
    end
for l = 1:2
    if plotting==1
    PL = eval(['plot' num2str(l)]);
    set(PL, 'NextPlot', 'replace')
    % PL.Position = Pos{g};
    histogram(PL, [h(ismember({h.lang}, language{l})).(estimator)]);
    refresh(WWW)
    set(PL, 'NextPlot', 'add')
    tit = title(PL, ['Direction of SI: ' language{l}]);
    tit.FontSize = 12;
    plot(PL,[cutoffOVERALL(1) cutoffOVERALL(1)], [0 2000], 'LineWidth', 2)
    plot(PL,[cutoffOVERALL(2) cutoffOVERALL(2)], [0 2000], 'LineWidth', 2)


    str = { ['ERPs based on epochs with *' estimator '* values']...
            ['less than ' num2str(q1) ' quantile and'       ]...
            ['greater than ' num2str(q2) ' quantile'        ]};
    axes(PL);
    o1 = text(0,1700,str);
    o1.FontSize = 12;
    o1.FontWeight = 'bold';
    xl = xlabel(PL, 'load measure');
    xl.FontSize = 12;
    yl = ylabel(PL, 'number of epochs with specific load');
    yl.FontSize = 12;
    end

ID1 = ismember({h.lang}, language{l});
ID2 = ismember({h.subj}, subjectsofinterest);
a = find(ismember([h.tmplab], 'low') & ID1 & ID2);
low = reshape([h(a).erp], length(Xtab.times), []);
low_mu = mean(low,2);
low_SE = std(low,0,2)/sqrt(size(low,2));

b = find(ismember([h.tmplab], 'high') & ID1 & ID2);
high = reshape([h(b).erp], length(Xtab.times), []);
high_mu = mean(high, 2);
high_SE = std(high,0,2)/sqrt(size(high,2));

m = find(ismember([h.tmplab], 'mid') & ID1 & ID2);
mid = reshape([h(m).erp], length(Xtab.times), []);
mid_mu = mean(mid, 2);
mid_SE = std(mid,0,2)/sqrt(size(mid,2));

    if plotting==1
    PL = eval(['plot' num2str(3+l)]);
    set(PL, 'NextPlot', 'replace');
    plot(PL, Xtab.times, gaussfilt(low_mu , kernel_len)) % plot ERPs convolved with a Gaussian (SD=6ms)
    set(PL, 'NextPlot', 'add');
    plot(PL, Xtab.times, gaussfilt(high_mu, kernel_len))
    plot(PL, Xtab.times, gaussfilt(mid_mu, kernel_len))
    
    Xtab.([win '_' language{l} '_' 'low']) = gaussfilt(low_mu , kernel_len);
    Xtab.([win '_' language{l} '_' 'mid']) = gaussfilt(mid_mu , kernel_len);
    Xtab.([win '_' language{l} '_' 'high']) = gaussfilt(high_mu , kernel_len);
    csvpath = ['/Users/RomanKoshkin/Documents/R/ERPs' estimator '.csv'];
    writetable(Xtab, csvpath, 'Delimiter', ',');
    ylabel (PL, '\muV', 'FontSize', 14)
    xlabel(PL, 'time, ms', 'FontSize', 14)
    title(PL, ['ERPs filtered with a Gaussian kernel (\sigma = 6 ms)  ' language{l}])
    PL.YMinorGrid = 'on'; PL.XMinorGrid = 'on';

    set(PL, 'ylim',[-2 2])
    leg = legend(PL, strcat('low (', num2str(length(a)),' trials)'),...
    strcat('high WM (', num2str(length(b)), ' trials)')); % , strcat('med WM (', num2str(length(m)), ' trials)'));
    leg.FontSize = 12;
    leg.Location = 'southeast';
if win=='mu_N1'
    S.Vertices = [LB_N1*1000 -2; LB_N1*1000 2; UB_N1*1000 2; UB_N1*1000 -2];
else
    S.Vertices = [LB_P1*1000 -2; LB_P1*1000 2; UB_P1*1000 2; UB_P1*1000 -2];
end
    S.Parent = PL; S.Faces = [1 2 3 4]; S.EdgeColor = 'none'; S.FaceAlpha = 0.25;
    patch(S)
    end

tA = [h(a).(win)];
tB = [h(b).(win)];
% [h,p] = ttest2(tA, tB , 'Vartype', 'unequal');
% X2 = 0;
for SubX = 1:length(subjectsofinterest)
    for WML = 1:length(WMloadsofinterest)
        idx =...
            ismember({h.subj}, subjectsofinterest{SubX}) &...
            ismember({h.lang}, language{l}) &...
            ismember([h.tmplab], WMloadsofinterest{WML});
        X2(SubX, WML) = mean([h(idx).(win)]);    
    end
end
[H,p] = ttest(X2(:,3), X2(:,1));                        %!!!!!!!

    if plotting==1
    axes(plot6)
    Y = [mean(X2(:,1)) mean(X2(:,2)) mean(X2(:,3))];
    X = [1 2 3];
    dodge = (rand(1,3)-0.5)./3;
    err = [std(X2(:,1))/sqrt(length(X2)) std(X2(:,2))/sqrt(length(X2)) std(X2(:,3))/sqrt(length(X2))];
    errorbar(X+dodge,Y, err, '-s','MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red', 'LineWidth', 1.5)
    set(plot6, 'NextPlot', 'add')
    plot6.XTick = [1 2 3] + dodge;
    plot6.XTickLabels = WMloadsofinterest;
    plot6.XTickLabels = WMloadsofinterest;
    plot6.YGrid = 'on';
    xlim([0 4]);
    ylim([-1.2 0]);
    legend(plot6, language)
    end

C(ii, jj, l) = p;
E(ii, jj, l) = mean(X2(:,1,1) - X2(:,3,1));

    if plotting==1
    pow = power_calc(...
        -0.5,...
        std(tA),...
        length(tA),...
        std(tB),...
        length(tB),...
        0.05);
    req_samp = samp_size_calc(abs(mean(tA)-mean(tB)), 0.8, 0.05, std([tA tB]));
    str = ['Pval = ' num2str(p) char(10) 'Power = ' , num2str(pow) ' ,' char(10) 'Req. N for 0.8 power = ' , num2str(req_samp)];

    Txt = text(-50,1.5,str);
    Txt.FontSize = 14;
    Txt.Parent = PL;
    end
end

    if plotting==1
    set(plot6, 'NextPlot', 'replace')  
    if ~exist('QQQ')
        QQQ = figure;
    end
        figure(QQQ);
    plotlines(h, QQQ, win)
    figure(WWW);
    end
% compute and display 2x2 ANOVA:
% % % % x = h(find(...
% % % %     ismember([h.tmplab], WMloadsofinterest) & ...
% % % %     ismember({h.subj}, subjectsofinterest)...
% % % %     ));
% % % % x(find(cellfun(@isempty,{x.tmplab}))) = [];
% % % % dir = {x.lang}';
% % % % load = [x.tmplab]'; % !!! iscellstr!!!!!!!!!!!
% % % % Y = [x.(win)];
% % % % A = subplot(2,3,6);
% % % % A.Position = Pos{6};
% % % % pos = get(subplot(2,3,6),'position');
% % % % delete(subplot(2,3,6))
% % % % [p, tbl] = anovan(Y, {dir, load},'model','full','varnames',{'direction','load'}, 'display', 'off');
% % % % t = uitable('Data', tbl(2:end, 2:end), 'ColumnName', tbl(1,2:end), 'RowName', tbl(2:end,1));
% % % % set(t,'units','normalized')
% % % % set(t,'position',Pos{6})
% % % % t.ColumnWidth = {64 40 30 55 55 50};

        if plotting==1
        str = [num2str(BL) ' of ' num2str(bigLoop) ' GRAND AVERAGE of ' num2str(length(subjectsofinterest)) ' subjects: ' strjoin(subjectsofinterest)];
        WWW.Name = str;
        set(WWW,'units','normalized','outerposition',[0 0 1 1]) % maximize figure

        drawnow
        set(QQQ,'WindowStyle','modal')
        end
    end
    rrr = toc;
    rrr
    disp([num2str(BL) ' of ' num2str(bigLoop) '|| ETA: ' num2str(rrr * (bigLoop-BL))])
end

% save data for max_projection_fork1.m:
% save('/Users/RomanKoshkin/Downloads/dataframe.mat', 'h')

% save data to R:
if saveR == 1
    h = rmfield(h, 'fullEp');
    h = rmfield(h, 'erp');
    writetable(struct2table(h), '/Users/RomanKoshkin/Documents/R/dataframe.csv', 'Delimiter', ',')
end
%%
figure
subplot(2,2,1)
A = imagesc(flipud(squeeze(C(:,:,1))));
ylabel('lower quantile');
xlabel('uppper quantile');
title(language{1});
A.Parent.YTick = [1:length(Q2)];
A.Parent.XTick = [1:length(Q2)];
A.Parent.YTickLabel = flip(Q2);
A.Parent.XTickLabel = Q2;
colorbar

subplot(2,2,2)
B = imagesc(flipud(squeeze(C(:,:,2))));
B.Parent.YTick = [1:length(Q2)];
B.Parent.XTick = [1:length(Q2)];
B.Parent.YTickLabel = flip(Q2);
B.Parent.XTickLabel = Q2;
title(language{2});
ylabel('lower quantile');
xlabel('upppr quantile');
colorbar

%%
subplot(2,2,3)
A = imagesc(flipud(squeeze(E(:,:,1))));
ylabel('lower quantile');
xlabel('uppper quantile');
title(['effect size' language{1}]);
A.Parent.YTick = [1:length(Q2)];
A.Parent.XTick = [1:length(Q2)];
A.Parent.YTickLabel = flip(Q2);
A.Parent.XTickLabel = Q2;
colorbar

subplot(2,2,4)
B = imagesc(flipud(squeeze(E(:,:,2))));
B.Parent.YTick = [1:length(Q2)];
B.Parent.XTick = [1:length(Q2)];
B.Parent.YTickLabel = flip(Q2);
B.Parent.XTickLabel = Q2;
title(['effect size' language{2}]);
ylabel('lower quantile');
xlabel('upppr quantile');
colorbar