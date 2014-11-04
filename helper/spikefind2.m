function [spikes] = spikefind2(x,thresh,dT)
% spikes = spikefind(x,thresh,dT)
%
% spikefind.m: find "spikes" in a given signal x
% x: input data
% thresh: amplitude threshold
% dT: minimum time above threshold (in # of time bins)
% spikes: indices (in x) where spikes occur
%
%
spikes=[];

ind = find(x>thresh);
di = diff(ind);
jumps = [0 find(di>1)];
if length(jumps)==1, return, end
tOn = ind(jumps+1); tOn=tOn(1:end-1);
tOff = ind(jumps(2:end));

for j = 1:length(tOn);
  if tOff(j)-tOn(j)>dT
    [~,j0] = max(x(tOn(j):tOff(j)));
    spikes = [spikes tOn(j)+j0-1];
  end
end

end