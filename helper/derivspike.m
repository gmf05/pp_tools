function [spikeInd] = derivspike(D,d2thresh)
% function [spikeInd,amp,derivs,secondDerivs] = derivspike(D,d2thresh)

% get normalized 1st, 2nd derivatives
dD = diff(D); %dD = dD./max(dD);
d2D = diff(dD); d2D = d2D./max(d2D);

% find - to + zero crossings of first deriv = local minima
minInd = find(dD(1:end-1)<0 & dD(2:end)>=0);

% keep places where 2nd deriv is > d2thresh
spikeInd = [];
for i = 1:length(minInd)
  if d2D(minInd(i))>d2thresh, spikeInd(end+1) = minInd(i); end;
end

% amps = D(spikeInd);
% derivs = dD(spikeInd);
% secondDerivs = d2D(spikeInd);

end