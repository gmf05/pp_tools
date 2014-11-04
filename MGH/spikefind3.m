function [spikes] = spikefind3(x,thresh1,dT,thresh2)
% spikes = spikefind(x,thresh1,dT,thresh2)
%
% spikefind.m: find "spikes" in a given signal x
% x: input data
% thresh1: amplitude threshold #1
% dT: minimum # of bins for which x > thresh1
% thresh2: minimum value of peak in x
%
% spikefind processes x in four steps:
% 1. find intervals where x > thresh1
% 2. keep intervals at least dTearticular intervals
% 4. keep times of peaks which are > thresh2
%
spikes=[];

ind = find(x>thresh1);
di = diff(ind);
jumps = [0 find(di>1)];
if length(jumps)==1, return, end
tOn = ind(jumps+1); tOn=tOn(1:end-1);
tOff = ind(jumps(2:end));

for j = 1:length(tOn);
  if tOff(j)-tOn(j)>dT
    [~,j0] = max(x(tOn(j):tOff(j)));
    j1 = tOn(j)+j0-1;
    if x(j1)>=thresh2
      spikes = [spikes tOn(j)+j0-1];
    end
  end
end

end