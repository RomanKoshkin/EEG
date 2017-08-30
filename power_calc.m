% this function computes observed power, and the minimum number of samples 
% to achieve the power desired:

function pow = power_calc(min_effect_size_wanted, sd1, n1, sd2, n2, alpha)
    % calculate the statistical power of a test
    % we want to find any effect that is larger than observed_effect_size
    % what is the power of the test that can detect this effect?
    SE = sqrt(sd1^2/n1 + sd2^2/n2);
    z = abs(norminv(alpha/2, 0, SE)); % compute the z_crit.
    Z = (z - abs(min_effect_size_wanted))/SE;
    pow = 1 - normcdf(Z, 0, 1);
end