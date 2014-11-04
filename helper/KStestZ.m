function [ks_stat,ks_ci,ks_p] = KStestZ(z)
% KStestZ.m
% [ks_stat, ks_ci] = ks_test(z)
% Input: z (rescaled ISIs), assuming exponential ISI rescaled -> uniform
% Output: ks_stat (KS score), ks_ci (confidence bound)
%
numISIs = length(z);

if numISIs<3
  error('Too few data points');
end

%---- rescale function is expontential               
% rs_fn = @(x)(1-exp(-x));
rs_cdf = @(x)(unifcdf(x,0,1));
% % %---- rescale function is identity
% % rs_fn = @(x)(x);
% % rs_cdf = @(x)(expcdf(x,1));

% z = rs_fn(z);
[eCDF,xCDF] = ecdf(sort(z));
aCDF = rs_cdf(xCDF);
ks_stat = max(abs(aCDF-eCDF));
ks_ci = 1.36/sqrt(numISIs+1);

n1     =  length(aCDF);
n2     =  length(eCDF);
n      =  n1 * n2 /(n1 + n2);
lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * ks_stat , 0);            
% 1-sided test:
% % %    ks_p  =  exp(-2 * lambda * lambda);
% 2-sided test (default):
%  "Use the asymptotic Q-function to approximate the 2-sided P-value."
j       =  (1:101)';
% ks_p  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2)); %?????
ks_p  =  1 - 2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2)); %???
ks_p  =  min(max(ks_p, 0), 1);


end