function KSplotZ(z,params)
% KStest.m
% [ks_stat, ks_ci] = ks_test(y,cif)
% Input: y (point process data), cif (cond'l intensity from model)
% Output: ks_stat (KS score), ks_ci (confidence bound)
%

if nargin<3, params.rs = 'exp'; params.alpha = 0.05; end

numISIs = length(z);

if numISIs<3
  error('Too few data points');
end

switch params.rs
  %---- rescale function is expontential               
  case 'exp'
    rs_fn = @(x)(1-exp(-x));
    rs_cdf = @(x)(unifcdf(x,0,1));
  %---- rescale function is identity
  case 'identity'
    rs_fn = @(x)(x);
    rs_cdf = @(x)(expcdf(x,1));
end

% z = rs_fn(z);
[eCDF,xCDF] = ecdf(sort(z));
aCDF = rs_cdf(xCDF);
% ks_stat = max(abs(aCDF-eCDF));
ks_ci = 1.36/sqrt(numISIs+1); % NOTE: get alpha level from params as well!

plot(eCDF,aCDF); hold on
x = [0:0.1:1];
plot(x,x,'r--',x,x+ks_ci,'r--',x,x-ks_ci,'r--','linewidth',1.5);

end