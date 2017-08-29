% load repeatedmeas
% rm = fitrm(between,'y1-y8 ~ Group*Gender'); % + Age + IQ','WithinDesign',within);
% rm = fitrm(between,'y1-y8 ~ Group*Gender + Age + IQ','WithinDesign',within);
% ranovatbl = ranova(rm)
% [ranovatbl] = ranova(rm,'WithinModel','w1+w2')
load ('myanova.mat')
for i = 1:1%length(between.subj)
    subject = char(between.subj(i));
    for j = 1:height(within)
        idx = find(...
            ismember({s.subj}, subject) &...
            ismember({s.lang}, char(within.w1(j))) &...
            ismember([s.tmplab], char(within.w2(j))));
        between(i,4+j) = table(mean([s(idx).mu_stat]));
    end
end

tab = table2array(between(:,5:10));
plot(mean(tab(:,1:3), 1))
hold on
plot(mean(tab(:,4:6), 1))

