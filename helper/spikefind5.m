function [spikes] = spikefind5(x,thresh1,dT1,thresh2,dT2,dT3)
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
% % % spikes=[];
% % % T = length(x);
% % % T0 = T-dT2;
% % % t=dT1+1;
% % % lastspike=-1e4;
% % % 
% % % % time loop:
% % % % march through time steps and check whether
% % % % left-hand-slope > thresh1 & right-hand slope < thresh2 
% % % % if so, local max is the spike
% % % while t<T0  
% % %   if x(t)-x(t-dT1)>=thresh1 && x(t)-x(t+dT2)>=thresh2 && t-lastspike>=dT3
% % %     % find local max:
% % %     [~,t0] = max(x(t-dT1:t+dT2));
% % %     t1 = t-dT1-1+t0;      
% % %     spikes = [spikes t1];
% % %     t=t+dT2;
% % %     lastspike = t1;
% % %   else
% % %     t=t+1; % no spike
% % %   end
% % % end

spikes=[];
T = length(x);
T0 = T-dT2;
t=dT1+1;
lastspike=-1e4;

% time loop:
% march through time steps and check whether
% left-hand-slope > thresh1 & right-hand slope < thresh2 
% if so, local max is the spike
while t<T0  
  if x(t)-x(t-dT1)>=thresh1 && x(t)-x(t+dT2)>=thresh2
    % find local max:
    [~,t0] = max(x(t-dT1:t+dT2));
    t1 = t-dT1-1+t0;      
      if t1-lastspike<=dT3
       if x(t1)>x(lastspike), lastspike=t1; spikes(end)=t1; end
       t=t+dT2;
      else
        spikes = [spikes t1];
        t=t+dT2;    
        lastspike = t1;
      end
  else
    t=t+1; % no spike
  end
end

end