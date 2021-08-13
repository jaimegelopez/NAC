function [corr12,mean1,mean2,var1,var2] = compute_continuous_variance(s1,s2,t)

%This script computes means and variances of continuous signals from
%hybrid simulation

deltat = t(end) - t(1);
mean1 = trapz(t,s1)/deltat;
mean2 = trapz(t,s2)/deltat;

s1_dev = s1 - mean1;
s2_dev = s2 - mean2;

var1 = trapz(t,s1_dev.^2)/deltat;
var2 = trapz(t,s2_dev.^2)/deltat;

cov12 = trapz(t,s1_dev.*s2_dev)/deltat;

corr12 = cov12/(sqrt(var1)*sqrt(var2));



end
