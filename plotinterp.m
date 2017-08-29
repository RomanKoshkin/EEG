function plotinterp(tempmt, subject, text_no)
est = [7 13]; % first is for non-interpolated, second for interpolated
xlpath = ['/Users/RomanKoshkin/Documents/MATLAB/EEG/' subject '_EEG/' subject '_arrows EEG_2.xls'];
xl_workbook = py.xlrd.open_workbook(xlpath);
xl_sheet = xl_workbook.sheet_by_index(int8(text_no-1));
qqq = cell(xl_sheet.col(int8(26)));
qqq = (cellfun(@char,qqq,'UniformOutput',false))';

figure
leng = 20;
str = qqq(1:leng);
for i = 1:length(str)
    str(i) = {strrep(str{i}, 'text:', '')};
    str(i) = {strrep(str{i}, '''', '')};
end
plot(tempmt(1:leng,2), tempmt(1:leng,est(1)));
ax1 = gca;
CL = tempmt(1:leng,2)+(0.001*(-1).^[1:leng])';
ax1.XTick = CL;
ax1.XTickLabel = str;
ax1.XTickLabelRotation = 90;
ax1.NextPlot = 'add';
ax1.YLim = [0 4];
ax1.YLabel.String = 'Weighted WM load';
ax1.FontSize = 14;
ax1.XTickLabelRotation = 45;
grid on
hold on

ax2 = axes('Position',ax1.Position, 'XAxisLocation','top', 'YAxisLocation','right', 'Color','none');
ax2.YTick = [];
ax2.YLim = [0 4];
ax2.XColor = 'r';
ax2.XGrid = 'on';
ax2.XLim = ax1.XLim;
plot(ax1, tempmt(find(tempmt(1:leng,8)),8), tempmt(find(tempmt(1:leng,est(2))),est(2)), 'LineWidth', 3, 'Color', 'r')
ax2.XTick = tempmt(find(tempmt(1:leng,8)),8);
ax2.XTickLabel = tempmt(find(tempmt(1:leng,8)),8);
ax2.XTickLabelRotation = 90;
ax2.GridLineStyle = '-.';
ax2.XLabel.String = 'Time, s';
ax2.FontSize = 14;
ax2.XTickLabelRotation = 45;
end
