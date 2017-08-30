% This functino calculates the require sample size to achieve the desired statistical
% power (https://www.coursera.org/learn/inferential-statistics-intro/lecture/kdnQf/power)
    
function n = samp_size_calc(d, desired_pow, desired_alpha, s)
    SE = d/(norminv(desired_pow,0,1)-norminv(desired_alpha/2,0,1));
    n = round(2*s^2/SE^2);
end
