function [ks_stat, ks_ci, z] = ks_test(y,cif)
% KStest.m
% [ks_stat, ks_ci] = ks_test(y,cif)
% Input: y (point process data), cif (cond'l intensity from model)
% Output: ks_stat (KS score), ks_ci (confidence bound)
%
spike_ind = find(y);           
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
ks_ci = 1.96/sqrt(numISIs+1);

end