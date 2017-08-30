% This script extracts relevant WM load data and corresponding timecodes
% from the Excel spreadsheets using the MATLAB-Python interface

%  WE TAKE DATA FROM THE CORRECTED CLRED _arrows EEG_3.xls:

subj = {'KOZ'}%{'ELT', 'SHE', 'BUL', 'POG', 'ROM', 'KOZ', 'KOS'}; % {'GRU'} % for GRU the sheet range is 0:3
for j = 1:length(subj)
    xlpath = ['/Users/RomanKoshkin/Documents/MATLAB/EEG/' subj{j} '_EEG/' subj{j} '_arrows EEG_3.xls'];
    EVSpath = ['/Users/RomanKoshkin/Documents/MATLAB/EEG/' subj{j} '_EEG/TimeCodeEVSdata_' subj{j} '.mat'];
    load(EVSpath)
    xl_workbook = py.xlrd.open_workbook(xlpath);
    tic
    for sheet_no = 0:7

        xl_sheet = xl_workbook.sheet_by_index(int8(sheet_no));
        column1 = cell(xl_sheet.col(int8(16))); %CWnored
        column2 = cell(xl_sheet.col(int8(28))); %timeIN
        column3 = cell(xl_sheet.col(int8(29))); %timeOFF
        column4 = cell(xl_sheet.col(int8(34))); %SYL
        column5 = cell(xl_sheet.col(int8(27))); %CWred
        column6 = cell(xl_sheet.col(int8(41))); %CLnored
        column7 = cell(xl_sheet.col(int8(42))); %CLred
        for i = 1:length(column1)
            if isfloat (column1{i}.value) & isfloat(column2{i}.value)
                values(i,1) = column1{i}.value;
                values(i,2) = column2{i}.value;
                values(i,3) = column3{i}.value;
                values(i,4) = column4{i}.value;
                values(i,5) = column5{i}.value;
                values(i,6) = column6{i}.value;
                values(i,7) = column7{i}.value;
            end
        end
        s.(['text' num2str(sheet_no+1)]) = values;
        values = 0;
    end
    save(EVSpath, 's')
    toc
end