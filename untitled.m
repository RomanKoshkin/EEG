% estimate = 'CLred';
% subjectsofinterest = {'POG', 'ELT', 'GRU', 'KOS', 'ROM', 'KOZ', 'SHE', 'BUL'};
% languageofinterest = {'er'};
% cutoff = [0.5 0.5];
% subjID = ismember({h.subj}, subjectsofinterest);
% langID = ismember({h.lang}, languageofinterest);
% IDofinteres = subjID & langID;
% Qlo = quantile([h(IDofinteres).(estimate)], cutoff(1));
% Qhi = quantile([h(IDofinteres).(estimate)], cutoff(2));
% a = find([h.(estimate)] < Qlo & IDofinteres); %quantile([h(IDofinteres).(estimate)], cutoff(1)));
% low = reshape([h(a).erp], 125, []);
% low_mu = mean(low,2);
% b = find([h.(estimate)] >= Qhi & IDofinteres); %quantile([h(IDofinteres).(estimate)], cutoff(2)));
% high = reshape([h(b).erp], 125, []);
% high_mu = mean(high,2);
% plot([-100:4:396],low_mu, 'LineWidth', 2);hold on; plot([-100:4:396],high_mu, 'LineWidth', 4);
% [quantile([h(IDofinteres).(estimate)], cutoff(1)) quantile([h(IDofinteres).(estimate)], cutoff(2))]
% [length(a), length(b)]