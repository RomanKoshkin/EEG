function ERPm = myplotX(epochs)
global EEG myERPplot ERPm 

epochs = find(~ismember(1:length(EEG.data), epochs));
channel = find(ismember({EEG.chanlocs.labels},'Cz') == 1);
y = squeeze (EEG.data(channel,:,epochs));
y = mean(y,2);
if length(y)==1
    y = squeeze (EEG.data(channel,:,epochs));
end
x = 1:125;
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