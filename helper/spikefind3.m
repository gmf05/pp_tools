function [spikes, slopesL, slopesR] = spikefind3(x,thresh,dT)
% spikes = spikefind(x,thresh1,dT,thresh2)
%
% spikefind.m: find "spikes" in a given signal x
% x: input data
% thresh1: amplitude threshold #1
% dT: minimum # of bins for which x > thresh1
% thresh2: minimum value of peak in x
%
% if deflection of at least thresh occurs over an interval dT
% then find local max following
%
if nargin<3
  spikeParams = thresh;
  thresh = spikeParams(1);
  dT = spikeParams(2);
end  

spikes=[];
slopesL=[];
slopesR=[];
T = length(x);
T0 = T-dT;
t=dT+1;

% time loop:
% march through time steps and check whether
% left-hand-slope > thresh1 & right-hand slope < thresh2 
% if so, local max is the spike
while t<T0  
  if x(t)-x(t-dT)>=thresh & x(t)-x(t+dT)>=thresh
    % find local max:
    [~,t0] = max(x(t-dT:t+dT));
    t1 = t-dT-1+t0;
    spikes = [spikes t1];
    if t1-dT<1
      mL = (x(t1) - x(1))/(t1-1);
    else
      mL = (x(t1) - x(t1-dT))/dT;
    end
    if t1>T0
      mR = (x(T) - x(t1))/(T-t1-1);
    else
      mR = (x(t1+dT) - x(t1))/dT;
    end
    slopesL = [slopesL mL];
    slopesR = [slopesR mR];
    t=t+dT;
  else
    t=t+1; % no spike
  end
end

end