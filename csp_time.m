% Data_c1 (c1) -- Ns x (T*Ntr)
% subtract the Data_c1 from Data_c2 (diff = Data_c2-Data_c1)
% dur -- length of time window in time points (not in ms)
% int -- step between time windows (not in ms)

function [dd1, finscore] =  csp_time(Data_c1,Data_c2,Fs,dur,int,start,stop,show)

[Ns T1] = size(Data_c1);
aC1 = mean(reshape(Data_c1, [Ns, Fs, T1/Fs]),3);
[Ns T2] = size(Data_c2);
aC2 = mean(reshape(Data_c2, [Ns, Fs, T2/Fs]),3);

lambda1 = 0;

clear dd1
clear finscore

h = 1;

% start>0, stop< Fs-dur, not in ms

for a = [start:int:stop]
    MMNRange = a:(a+dur);
    OtherRange = [1:(a-1),(a+dur+1):500];

    DiffStdDev1 = -aC1+aC2;
    C1 = DiffStdDev1(:,MMNRange)*DiffStdDev1(:,MMNRange)'/size(DiffStdDev1,1);
    C2 = DiffStdDev1(:,OtherRange)*DiffStdDev1(:,OtherRange)'/size(DiffStdDev1,1);

    [v1,d1] = eig(C1+0.05*trace(C1)/Ns*eye(size(C1)),C2+0.05*trace(C2)/Ns*eye(size(C2)));
    [~, order] = sort(diag(d1),'descend');
    v1 = v1(:,order);
    v1 = v1*diag(1./sqrt(sum(v1.^2,1)));

    % проверить знак
    prStd1 = v1(:,1)'*aC1;
    prDev1 = v1(:,1)'*aC2;
    [m ind] = max(abs(prDev1(MMNRange)));
    if prDev1(MMNRange(ind))<0
        v1(:,1) = -v1(:,1);
    end

    dd1(h)= max(max(d1));
    
    
    % permutation test
   
    AllD1 = [Data_c1, Data_c2];
    Ntr1 = size(AllD1,2)/Fs;

    Niter = 200;
    
    ss_true(1,:) = v1(:,1)'*(-aC1+aC2);
    true_score = abs(mean(ss_true(MMNRange))-mean(ss_true(OtherRange)));
    
    for j = 1:Niter
        ind1 = randperm(Ntr1);
        as_perm1 = zeros(Ns, Fs);
        ad_perm1 = zeros(Ns, Fs);
        
        for i = 1:(T1/Fs)
            as_perm1 = as_perm1 + AllD1(:, (Fs*(ind1(i)-1)+1):ind1(i)*Fs);
        end
        for i = (T1/Fs+1):Ntr1
            ad_perm1 = ad_perm1 + AllD1(:, (Fs*(ind1(i)-1)+1):ind1(i)*Fs);
        end
        as_perm1 = as_perm1/(T1/Fs);
        ad_perm1 = ad_perm1/(Ntr1-T1/Fs);

        DiffStdDev1 = -as_perm1+ad_perm1;
        C11 = DiffStdDev1(:,MMNRange)*DiffStdDev1(:,MMNRange)'/size(DiffStdDev1,1);
        C21 = DiffStdDev1(:,OtherRange)*DiffStdDev1(:,OtherRange)'/size(DiffStdDev1,1);

        [v_perm1,d_perm1] = eig(C11+0.05*trace(C11)/Ns*eye(size(C11))-lambda1*eye(size(C11,1), size(C11,2)),C21+0.05*trace(C21)/Ns*eye(size(C21)));
        [~, order] = sort(diag(d_perm1),'descend');
        v_perm1 = v_perm1(:,order);
        d_perm1 = d_perm1(:,order);
        v_perm1 = v_perm1*diag(1./sqrt(sum(v_perm1.^2,1)));


        ss_hat(1,:)= v_perm1(:,1)'*(-as_perm1+ad_perm1);
        score(j) = abs(mean(ss_hat(MMNRange))-mean(ss_hat(OtherRange)));
    end
       
    finscore(h) = sum(score>true_score)/Niter;
    
    if show == 1
        [finscore(h), a, h]
    end
    h = h+1;
end       

end