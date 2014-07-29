function [spikeInd,dD,d2D] = derivspike(D,d1thresh,d2thresh)
% function [spikeInd,amp,derivs,secondDerivs] = derivspike(D,d1thresh,d2thresh)

% get 1st & 2nd derivatives
dD = zscore(diff(D));
d2D = zscore(diff(dD));

% find - to + zero crossings of first deriv = local minima
minInd = find(dD(1:end-1)<0 & dD(2:end)>=0);

% keep places where 2nd deriv is > d2thresh
spikeInd = [];
for i = 2:length(minInd)-1
  if min(dD(minInd(i-1):minInd(i)))<-d1thresh && max(dD(minInd(i):minInd(i+1)))>d1thresh && ...
      d2D(minInd(i))>d2thresh
    spikeInd(end+1) = minInd(i);
  end
end

end