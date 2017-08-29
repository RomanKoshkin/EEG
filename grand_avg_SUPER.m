function grand_avg_SUPER(subject, est_vec, LB, UB)
load('params.mat')
clear AccC1 AccC2

times = [-100:4:396];
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
%             s(j).subj = subject{i};
%             s(j).lang = language{i};
%             s(j).WMloadEst = estimate{i};
%             s(j).erp = 
        end
            
            C1 = mean(AccC1,2);
            C2 = mean(AccC2,2);

             
            subplot (length(language),length(estimate), u)
            plot(times, C1, times, C2)
            grid on
            q = title(['Units: ' estimate{k} ', Language: ' language{l}]);
            q.FontSize = 13;
            u = u+1;
            legend('LOW WM', 'high WM')
            
            
            M1 = mean(C1);
            M2 = mean(C2);
            
            [h,p] = ttest2(C1(LB:UB),C2(LB:UB), 'Vartype', 'unequal');
            
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
            

            S.Vertices = [Lb*1000 -2; Lb*1000 2; Ub*1000 2; Ub*1000 -2];
            S.Faces = [1 2 3 4];
            S.EdgeColor = 'none';
            S.FaceAlpha = 0.25;
            patch(S)
            
            AccC1 = []; AccC2 = [];
    end
end
            str = ['GRAND AVERAGE of ' num2str(length(subject)) ' subjects: ' strjoin(subject)];
            ylim=get(gcf,'ylim');
            xlim=get(gcf,'xlim');
            q1 = text(xlim(end)/2,1.5,str);
            q1.FontSize = 16;
            q1.FontWeight = 'bold';
            q1.Position = [-322 14 0];
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) % maximize figure
end