clf;
newtimef(low, EEG.pnts, [-100 396], EEG.srate,... %ALLEEG(2).data(22,:,:)
    'cycles', [1 0.5],...   % number of cycles (of the wavelet) per window
    'winsize',50,...        % small windows will not capture long waves
    'padratio', 4,...       % kind of increases resolution in the frequency domain
    'maxfreq', 45);%,...       % maximum frequency plotted
%     'alpha', 0.05,...       $ permutation test parameters
%     'naccu', 200);         


% sampling rate = 250 Hz
% maximum frequency captured by this sRate = 250/2 % by the Nyquist Theorem
% suppose we have a wave with freq 10 Hz, how many samples is a full cycle?
% 
% 1 Hz will need a 250-sample window,
% 10 Hz will need a 25 sample window