function [ks_stat, ks_ci, z, ks_p] = KStest(y,cif)
% KStest.m
% [ks_stat, ks_ci] = ks_test(y,cif)
% Input: y (point process data), cif (cond'l intensity from model)
% Output: ks_stat (KS score), ks_ci (confidence bound)
%

spike_ind = [0; find(y)];
numISIs = length(spike_ind)-1;

if numISIs<3
  error('Too few data points');
end

z = zeros(1, numISIs);
for j=1:numISIs                           
  z(j) = sum(cif(spike_ind(j)+1:spike_ind(j+1)));
end

%---- rescale function is expontential               
rs_fn = @(x)(1-exp(-x));
rs_cdf = @(x)(unifcdf(x,0,1));
% % %---- rescale function is identity
% % rs_fn = @(x)(x);
% % rs_cdf = @(x)(expcdf(x,1));

z = rs_fn(z);
[eCDF,xCDF] = ecdf(sort(z));
aCDF = rs_cdf(xCDF);
ks_stat = max(abs(aCDF-eCDF));
ks_ci = 1.36/sqrt(numISIs);

% % % ks_p = NaN;
% compute ks_p: p-value for ks-statistic
% WARNING: still under construction
n1     =  length(aCDF);
n2     =  length(eCDF);
n      =  n1 * n2 /(n1 + n2);

lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * ks_stat , 0);            
% 1-sided test:
% % %    ks_p  =  exp(-2 * lambda * lambda);
% 2-sided test (default):
%  "Use the asymptotic Q-function to approximate the 2-sided P-value."
j       =  (1:101)';
% ks_p  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
ks_p  =  1 - 2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
ks_p  =  min(max(ks_p, 0), 1);

%
% cinv = @(x)(x);
% ks_p = cinv(ks_stat*sqrt(n));

end