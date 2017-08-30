% Run this script first thing on the raw dataset. Renames event codes 
% (set automatically by the BrainVision amplifier) to simple numeric values
% to match the text numbers.

j = 1;
for q = 1:8 % specify the number of text event lables (e.g. from 'S  1':'S  8')
            % i.e. in my experiment there was one single EEG file
            % containing EIGHT conditions. Each condition had its unique
            % trigger (S  1, S  2 etc.)
    while j <= length(EEG.event)
        if strcmp(EEG.event(j).type, ['S  ',num2str(q)])
            EEG.event(j).type = num2str(q);
        end
        j = j + 1;
    end
    j = 1;
end
clear q j