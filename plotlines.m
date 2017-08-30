% This function creates subject-specific interaction plots. Normally called
% from the NEW_GRAND_AVERAGE.m

function plotlines (s, fig, win)
    subj = {'KOK', 'ELT', 'KOZ', 'POG', 'KOS', 'ROM', 'SHE', 'BUL','GRU'};
    figure(fig);
    if strcmp(get(fig,'CurrentCharacter'), 'q')==1
        tmp = input('User Break');
    end

    idx1 = ismember({s.lang}, 'er');
    idx2 = ismember({s.lang}, 'RE');
    idx3 = ismember([s.tmplab], 'high');
    idx4 = ismember([s.tmplab], 'mid');
    idx5 = ismember([s.tmplab], 'low');
    X = [1 2 3];

    for l = 1:length(subj)
        dodge = (rand(1,3)-0.5)./15;
        idxSub = ismember({s.subj}, subj{l});
        E1 = [s(idxSub & idx1 & idx3).(win)];
        E1mu = mean(E1);
        E1se = std(E1)/sqrt(length(E1));
            R1 = [s(idxSub & idx2 & idx3).(win)];
            R1mu = mean(R1);
            R1se = std(R1)/sqrt(length(R1));

        E2 = [s(idxSub & idx1 & idx4).(win)];
        E2mu = mean(E2);
        E2se = std(E2)/sqrt(length(E2));
            R2 = [s(idxSub & idx2 & idx4).(win)];
            R2mu = mean(R2);
            R2se = std(R2)/sqrt(length(R2));

        E3 = [s(idxSub & idx1 & idx5).(win)];
        E3mu = mean(E3);
        E3se = std(E3)/sqrt(length(E3));
            R3 = [s(idxSub & idx2 & idx5).(win)];
            R3mu = mean(R3);
            R3se = std(R3)/sqrt(length(R3));

        if isnan(E2mu) == 1 || isnan(R2mu) == 1
            E2 = 0;
            R2 = 0;
        end

        ww = subplot (3, 4, l);
        errorbar(X+dodge, [E1mu E2mu E3mu], [E1se E2se E3se]); hold on
        dodge = (rand(1,3)-0.5)./8;
        errorbar(X+dodge, [R1mu R2mu R3mu], [R1se R2se R3se]); hold off
        xlim([0 4]);
        ylim([-2 1]);
        ww.XTick = X;
        ww.XTickLabels = {'high', 'mid', 'low'};
        legend({'er', 'RE'})
        title(subj{l})
    end
end
