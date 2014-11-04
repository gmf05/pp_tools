function y = cumdownsample(x,dsFactor)
N = length(1:dsFactor:size(x,2));
y = 0*x(:,1:N);
for n = 1:N-1
  y(:,n) = sum(x(:,(n-1)*dsFactor+(1:dsFactor)),2);
end
end