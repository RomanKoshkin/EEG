% This function takes epoch numbers as the only parameter. It finds these
% specified epochs in the current EEGlab dataset,
% averages them and plots the average ERP waveform. 
% The function is convenient for quickly exploring the different combination
% of epochs. Speeds up time compared to using EEGlab's GUI.

function ERPm = myplotX(epochs)
global EEG myERPplot ERPm 

epochs = find(~ismember(1:length(EEG.data), epochs));
channel = find(ismember({EEG.chanlocs.labels},'Cz') == 1);
y = squeeze (EEG.data(channel,:,epochs));
y = mean(y,2);
if length(y)==1
    y = squeeze (EEG.data(channel,:,epochs));
end
x = 1:150;
try 
    figure(ERPm)
catch
    ERPm = figure;
end
figure(ERPm)
myERPplot = plot(x, y);
grid on
set(myERPplot,'XData',x,'YData',y);
end