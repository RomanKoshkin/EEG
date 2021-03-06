% This is a legacy function that plots grand average ERPs by subject. We
% used the improved NEW_GRAND_AVERAGE.m instead.

function grand_avg_SUPER(subject, est_vec, LB, UB)

    % load boundary values for different WM load estimators:
    load('params.mat')

    clear AccC1 AccC2

    % initialize times, sampling rate
    times = [-100:4:496];
    xmin = -0.1;
    srate = 250;

    Lb = LB;
    Ub = UB;
    LB = round(-xmin*srate+LB*srate);
    UB = round(-xmin*srate+UB*srate);
    load ('AV.mat')

    AccC1 = [];
    AccC2 = [];
    k = 0;
    u = 1;

    language = {'both', 'ER', 're'};
    estimate = EST(est_vec);
    for l = 1:length(language)
        for k = 1:length(estimate) 
            for i = 1:length(subject)
                AccC1 = [AccC1 AV.(subject{i}).(language{l}).(estimate{k}).C1.data];
                AccC2 = [AccC2 AV.(subject{i}).(language{l}).(estimate{k}).C2.data];
            end

            % compute mean voltages:
            C1 = mean(AccC1,2);
            C2 = mean(AccC2,2);

            % plot grand average ERPs comparing the low and high WM load
            % conditions defined using fixed cutoffs (params.mat):
            subplot (length(language),length(estimate), u)
            plot(times, C1, times, C2)
            grid on
            q = title(['Units: ' estimate{k} ', Language: ' language{l}]);
            q.FontSize = 13;
            u = u+1;
            legend('LOW WM', 'high WM')


            M1 = mean(C1);
            M2 = mean(C2);

            % perform an independent samples t-test (exploratory).
            [h,p] = ttest2(C1(LB:UB),C2(LB:UB), 'Vartype', 'unequal');

            % add legends, captions and additional info to the plot:
            str = strcat('MeanLow=', num2str(M1),', MeanHigh=', num2str(M2),char(10), ', Pval=', num2str(p));
            ylim=get(gca,'ylim');
            xlim=get(gca,'xlim');
            q = text(xlim(1)+5,1.8,str);
            q.FontSize = 12;
            str = ['trials in LO=' num2str(length(AccC1)) ', trials in HI=' num2str(length(AccC2))];
            ylim=get(gca,'ylim');
            xlim=get(gca,'xlim');
            q = text(xlim(1)+5,1.5,str);
            q.FontSize = 12;

            if k == 1 
                str = ['WM estimator: ' num2str(22) '. (1 - compression, 5 - no compression)'];
                ylim=get(gca,'ylim');
                xlim=get(gca,'xlim');
                q = text(xlim(1)+5,1.2,str);
                q.FontSize = 12;
            end

            % highlight the window of interest:
            S.Vertices = [Lb*1000 -2; Lb*1000 2; Ub*1000 2; Ub*1000 -2];
            S.Faces = [1 2 3 4];
            S.EdgeColor = 'none';
            S.FaceAlpha = 0.25;
            patch(S)

            AccC1 = []; AccC2 = [];
        end
    end
    
    % add a title for the entire figure with more useful info:
    str = ['GRAND AVERAGE of ' num2str(length(subject)) ' subjects: ' strjoin(subject)];
    ylim=get(gcf,'ylim');
    xlim=get(gcf,'xlim');
    q1 = text(xlim(end)/2,1.5,str);
    q1.FontSize = 16;
    q1.FontWeight = 'bold';
    q1.Position = [-322 14 0];
    set(gcf,'units','normalized','outerposition',[0 0 1 1]) % maximize figure
end