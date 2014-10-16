function spikes = spikefind(x,thresh,dT)
% spikes = spikefind(x,thresh,dT)
%
% spikefind.m: find "spikes" (steep extrema) of a given signal x
% x: input data
% thresh: "steepness" threshold (where steepness = product of derivatives)
% dT: time window on each side (number of time bins)
% spikes: indices (in x) where spikes occur
%
% sample values, e.g.
% dT = round(0.03*dt);
% thresh = -0.04;

% might want to zscore x & adjust thresh, if x has large amplitudes
% if range(x)>1e2, x = zscore(x); end

T = length(x);
spikes=[];

for t = dT+1:dT:T-dT
  t1 = t-dT;
  t2 = t+dT;
  dx1=x(t)-x(t1);
  dx2=x(t2)-x(t);
  dx1*dx2
  if dx1*dx2<thresh && dx1<dx2    
    [~,ti] = min(x(t1:t2)); % find local min on [t1,t2]
    t0 = t1-1+ti; % index shifted from [t1,t2] to [1,T]
    spikes=[spikes t0];
  end
end

end