subject = 'KOZ';
sheet_no = int8(0);
path_to_file = ['/Users/RomanKoshkin/Documents/MATLAB/EEG/' subject '_EEG/']
file = [subject '_arrows EEG_2.xls'];
xl_workbook = py.xlrd.open_workbook([path_to_file file]);
xl_sheet = xl_workbook.sheet_by_index(sheet_no);
strsplit(char(xl_sheet.name), '#')
% if you need to export the structure as a csv file:
% struct2csv(s,'test.csv')

for i = 1:xl_sheet.nrows
    b{i} = char(column1{i}.value);
end
data = [EEG.event.latency]./EEG.srate;
edges = [b{:,2}; b{end,3}]';
Q = discretize(data, edges);
B = cell2struct(b,{'q', 'qwe', 'qwer'},2);
tic
for i = 1:length(edges)
    b{i,4} = [EEG.event(find(Q==i)).latency]/EEG.srate;
    b{i,5} = ((b{i,4} - b{i,2})/(b{i,3} - b{i,2}) - 0.5) * 2;
end
toc
    