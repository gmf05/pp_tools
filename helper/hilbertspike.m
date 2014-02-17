function [spike_ind] = hilbertspike(d0, threshold, min_refract)
%
%   [spike_times] = hilbertspike2(d0, threshold, min_refract)
%
%   Written by Grant Fiddyment, Boston University, October 2011
%   
%   Summary: Finds "spikes" in the voltage trace d0 separated by at least
%   `threshold` by identifying the spikes in the Hilbert transform of the
%   signal.
%   
%    Description of inputs:
%    ------------------------------------------------------------
%           d0 = Voltage trace of a single electrode
%    threshold = A thresholding parameter describing what change in voltage
%                is deemed a spike. Namely, a we declare a spike if the
%                *peak-to-peak* change in voltage is at least `threshold`
%  min_refract = The minimum number of time steps allowed between spikes
%                Originally used in an attempt to eliminate line noise;
%                now mostly we choose min_refract = 1 (default)
%           dt = Time binning resolution
% [tmin, tmax] = Time window over which to find spikes
%
%   Description of outputs:
%   -------------------------------------------------------------
%   spike_times = A list of times when spikes occur
%   

% % threshold = 1*std(d0)

% if nargin<2, threshold = 100; end
if nargin <3, min_refract = 1; end
if nargin<4, dt=1; end
t_axis = (1:length(d0))*dt;
if nargin<5, tmin=0; end
if nargin<6; tmax=t_axis(end); end

window_span = round(4e-3/dt);

h_transform=imag(hilbert(d0));      % Hilbert transform of the signal
t0 = getclosest(t_axis,tmin);            % Index of the beginning time
t1 = getclosest(t_axis,tmax);            % Index of the end time
h_spikes=[];                        % Empty array for time indices of spikes

% First we find a sequence of local max/min pairs, max_ind and min_ind
htderiv = diff(h_transform(t0:t1));
temp1 = htderiv(1:end-1);
temp2 = htderiv(2:end);
ind = find(temp1.*temp2 <= 0);

max_ind = [];  min_ind = [];

for k=1:length(ind)
    if htderiv(ind(k)) > 0
        max_ind = [max_ind ind(k)+t0];
    else
        min_ind = [min_ind ind(k)+t0];
    end
end

N = min(length(max_ind), length(min_ind));
heights=[]; k=1;
count=1;

while k < min(length(max_ind), length(min_ind))

% Adjust the sequence of pairs so that it alternates between min/max
% (i.e. no max/max or min/min pairs) 
    while max_ind(k) > min_ind(k)        
        min_ind(k)=[];
    end    
    
    % Check whether the max/min pair represents a voltage change of at 
    % least `threshold`
    if h_transform(max_ind(k))-h_transform(min_ind(k)) >= threshold
            
        % Two options for identifying the spike between max/min:
        % (1) Take time halfway between time of max and time of min
        temp_ind = round((max_ind(k)+min_ind(k))/2);
        
        % (2) Take time when voltage is halfway between max and min
%         height = (h_transform(max_ind(k))+h_transform(min_ind(k)))/2;
%         temp_ind = getclosest(h_transform(max_ind(k):min_ind(k)), height)+max_ind(k);
        
%         window_span = 4; % how much "jitter" to allow b/w HT spikes and raw signal spikes
		ws = min(temp_ind-1, window_span);
        [~, temp_ind2] = min(d0(temp_ind-ws:temp_ind+ws));
        temp_ind = temp_ind-ws-1+temp_ind2;
                
        %If this spike is > min_refract away from last spike (OR this is the
        %first spike), save spike time
        if isempty(h_spikes) || temp_ind - h_spikes(end) >= min_refract
            h_spikes = [h_spikes temp_ind];
            count=count+1;
        end
    end
    k=k+1;
end

% Plot Hilbert transform spikes
%plot(t(h_spikes), d0(h_spikes), 'rx', 'MarkerSize', 8, 'LineWidth', 2);

%   Remove any duplicates
spike_ind = unique(h_spikes);

end


function inds=getclosest(data,vals)
        inds=zeros(1,length(vals));
        
        for i=1:length(vals)
            [y ind]=min(abs(data-vals(i)));
            inds(i)=ind;
        end

end
