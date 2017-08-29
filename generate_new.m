function [Data_s, Data_d, Ns, T] = generate_new(bGenNewNoise, Ntr, lambda, jit)
    
    % forward model matrix 
    G = load('FM2.mat');
    Gain = G.FM2.Gain;
    Gain = Gain([1:4, 6:24, 26:32, 34:62],:); % Gain = Gain([1:4, 6:24, 26:30],:); %
    Gx = Gain(:,1:3:size(Gain, 2));
    Gy = Gain(:,2:3:size(Gain, 2));
    Gz = Gain(:,3:3:size(Gain, 2));
    R = G.FM2.GridLoc;

    % location of four dipoles (on the both sides one for N1 and one for MMN)
    XYZGen = [0.005,0.06,0.05; 0.005, -0.06, 0.05; 0.01, 0.06, 0.055; 0.01, -0.06, 0.055]; % 0.01,-0.05,0.05; 0.015, 0.05, 0.055;  0.015, -0.05, 0.055; -0.05, 0, 0.07];
    % Ns --  number of sensors, Nsites -- number of sourses
    [Ns, Nsites] = size(Gx);

    % Find corresponding topographies
    % ??? ??????? ?????????
    for i=1:size(XYZGen,1)
        d = repmat(XYZGen(i,:),size(R,1),1) - R;
        % ?????????? ?? ????????? ?? ?????? ?? ??????
        d = sum(d.*d,2);
        % ??????? ????????? ??????? ? ??????????? ??????????
        [val ind] = min(d);
        % ?????????? ????????? ???????
        XYZGenAct(i,:) = R(ind,:);
        % take the first dipole in the tangent plane
        GenInd(i) = ind; 
    end;

    % data generation
    T = 500; % number of time points
    Fs = 500;

    % ?????????? ????????, ?????????? ????????? ????? ?? ???? ???????? 
    ERPd = zeros(Ns,Ntr*T);
    ERPs = zeros(Ns,Ntr*T);
    BrainNoise  = zeros(Ns,Ntr*T);
    SensorNoise  = zeros(Ns,Ntr*T);

    % ????????? ??? ??????? ????????? ??????????
    t = linspace(0,1,Fs); % ?????????????? ?????????? ?? 0 ?? 1, ? ????
    alpha = 0.25; % coupling strength reciprocal [0:1]
    F1 = 3; % Hz 
    F2 = 2; % Hz
    phi1= 0;

    Gz = Gz/100;
    Gx = Gx/100;
    Gy = Gy/100;

    
    range = 1:2;
    % for each node
    for i=1:Nsites
        % ?????? ? ??????? g (?????????? ???????? ? 3 ?????????? ?? ?,?,z) 
        g = [Gx(:,i) Gy(:,i) Gz(:,i)];
        % ??????? svd ?????????? ??????? g
        [u s v] = svd(g);
        G2d(:,range) = u(:,1:3)*v(:,[1,2]);
        range = range + 2;
    end;

    
    clear s
    range = 1:T;
    for tr = 1:Ntr
        jitter = randn()/jit*0.1;
        s(1,:) = sin(2*pi*F1*(t-jitter)+phi1).*exp(-10*(t-0.1-jitter).^2);
        s(2,:) = sin(2*pi*F1*(t-jitter)+phi1).*exp(-10*(t-0.1-jitter).^2);
        s(3,:) = -0.15*exp(-200*(t-0.25-jitter).^2);
        s(4,:) = -0.15*exp(-200*(t-0.25-jitter).^2); % MMN ????? 100-150

        Ggen1 = 0*Gx(:,GenInd(1))+ 0*Gy(:,GenInd(1)) + 1*Gz(:,GenInd(1));
        Ggen2 = 0*Gx(:,GenInd(2))+ 0*Gy(:,GenInd(2)) + 1*Gz(:,GenInd(2));
        Ggen3 = 0*Gx(:,GenInd(3))+ 0*Gy(:,GenInd(3)) + 1*Gz(:,GenInd(3));
        Ggen4 = 0*Gx(:,GenInd(4))+ 0*Gy(:,GenInd(4)) + 1*Gz(:,GenInd(4));

        % ??????-?? ????????? ??????? ????????? ???????? ?????? ?? ??????????
        % ?? ?, ? ??????? ?? ??????????
        erpd = ([Ggen1 Ggen2]*s(1:2,:)+[Ggen3 Ggen4]*s(3:4,:)); % Nsensors x T
        ERPd(:,range) = erpd;

        % ????????
        erps = [Ggen1 Ggen2]*s(1:2,:);
        ERPs(:,range) = erps;

        if (bGenNewNoise)
                brainnoise = GenerateBrainNoise(G2d,T,500,1000,Fs);
                brainnoise = brainnoise/sqrt(sum((brainnoise(:).^2)));
                BrainNoise_d(:,range) = brainnoise;
                brainnoise = GenerateBrainNoise(G2d,T,500,1000,Fs);
                brainnoise = brainnoise/sqrt(sum((brainnoise(:).^2)));
                BrainNoise_s(:,range) = brainnoise;
        end;
        sensornoise = randn(size(erpd));
        sensornoise = sensornoise/sqrt(sum((sensornoise(:).^2)));
        SensorNoise_s(:,range) = sensornoise;
        sensornoise = randn(size(erpd));
        sensornoise = sensornoise/sqrt(sum((sensornoise(:).^2)));
        SensorNoise_d(:,range) = sensornoise;
        range = range+T;
        %tr
    end

    if(~bGenNewNoise)
        load BrainNoiseFile_1000
    else
        save BrainNoiseFile_1000 BrainNoise_s BrainNoise_d;
    end

    Nnoise = size(BrainNoise_s,2)/Fs;
    
    ind_S = randi([1 (Nnoise-Ntr)],1,1);
    ind_D = randi([1 (Nnoise-Ntr)],1,1);
   
    BrainNoise_S = BrainNoise_s(:,((ind_S-1)*T+1):((ind_S-1)*T+T*Ntr));
    BrainNoise_D = BrainNoise_d(:,((ind_D-1)*T+1):((ind_D-1)*T+T*Ntr));
    
    clear Data_s Data_d;

    Data_s = lambda*BrainNoise_S + ERPs + 0*SensorNoise_s;
    Data_d = lambda*BrainNoise_D + ERPd + 0*SensorNoise_d;
    
end