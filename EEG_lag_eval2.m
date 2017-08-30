%% RUN THIS SCRIPT SECOND ON RAW DATA TO CORRECT THE EVENT TIMING
% THE ONLY INPUT IS THE k VARIABLE (NUMBER OF TEXTS IN THE EEG). EACH
% TEXT SHOULD HAVE ITS UNIQUE NUMERIC TRIGGER CODE FOR THE PROBES.
    
A = 0;
C = 0;

for k = 1:8
    K = num2str(k);
    j = 0; % initialize j to zero

    for i =  1:length(EEG.event)
        if strcmp(EEG.event(i).type, K) == 1    % EVENT OF INTEREST
            j = j + 1;
            A(j,1) = EEG.event(i).latency;      % extract trigger latencies (as measured in samples/s)
        end
    end
    
    A(2:end,2) = diff(A);                       % compute timing stability of event latencies
    C = abs(EEG.data(33,1:size(EEG.data,2)))>8000;        % extract samples where the probes sound was on (square voltage plateaus)
    
    %##############################################################################
    % !!!!!!!!!!!!!!!for FEB 04 EEG rec 6500 threshold for S'  1', 4000 for 'S  2'
    %##############################################################################
    j = 1;
    i = A(1,1);
    while i < length(C)
        if C(1,i) == 0 & C(1,i+1) == 1 & j <= size(A,1) % find real onsets of sound probes (in samples)

            A (j,3) = i;             % put those onsets into the comparison matrix
            if j < size(A,1)
                i = A(j+1,1);
            end
            j = j+1;
        end

        i = i + 1;
    end

    A(:,4) = A(:,3)-A(:,1);             % compute lags between trigger and probe onsets

    disp({'mean lag btwin triggr &' 'probe onset' mean(A(3:(size(A,1)-1),4))/EEG.srate 'seconds'})
    disp({'standard deviation' std(A(3:(size(A,1)-1),4))/EEG.srate 'seconds'})
    clear j i 

    for i = 1:length(A)
        for j = 1:length(EEG.event)
            if EEG.event(j).latency == A(i,1)
                EEG.event(j).latency = A(i,3);
            end
        end
    end


    disp ([K, ' Data corrected for lags between trigger and sound onsets'])
    clear A C
end
clear k i j K


% %% NOW FIND OFFSETS AND INSERT EVENT MARKERS INTO THE DATASET
% j=1;
% i=1;
% 
% for i = 1:length(A)
% E = 250;                           % expected duration of the probes (in samples)
%    while E > 1
%         if C(A(i,3) + E) == 0 & C(A(i,3) + E - 1) == 1 % find real offsets of sound probes (in samples)
%             A (i,5) = A(i,3) + E;   % put those offsets in column 5 the comparison matrix
%             E = 1;                  % exit while
%         end
%         E = E - 1;
%    end
% end
% 
% A(:,6) = A(:,5)-A(:,1);             % compute probe length
% 
% disp({'probe mean length' mean(A(3:(size(A,1)-1),6))/EEG.srate 'seconds'})
% disp({'probe SD' std(A(3:(size(A,1)-1),6))/EEG.srate 'seconds'})
% clear j i 
% 
% % now update the EEG.event structure
% k = length(EEG.event);
% for i = 1:size(A,1)-1
%     EEG.event(i+k).latency = A(i,5);
%     EEG.event(i+k).duration = 1;
%     EEG.event(i+k).channel = 0;
%     EEG.event(i+k).bvmknum = i+k;
%     EEG.event(i+k).type = 'S 3';
%     EEG.event(i+k).code = 'Stimulus';
%     EEG.event(i+k).urevent = i+k;
% end
% 
% disp ('Real offset times computed.')
%     
    