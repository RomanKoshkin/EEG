% Gaussian filter for ERP waveforms.

function filtered = gaussfilt (vec, sigmaMS)    % the function assumes sigma in ms
    samplelength = 4;                           % sample length, ms
    sigma = sigmaMS/samplelength;               % and converts it to samples
    mu = 0;
    sz = 125;                                   % length of gaussFilter vector
    x = linspace(-sz / 2, sz / 2, sz);
    gaussian = 1/(sigma*sqrt(2*pi))*exp(-(x-mu).^2/(2*sigma^2));
    gaussian = gaussian / sum (gaussian);       % normalize

    % plot (x, gaussian)
    filtered = conv (vec, gaussian, 'same');    % we convolve our vector with the
                                                % Gaussian kernel
end
