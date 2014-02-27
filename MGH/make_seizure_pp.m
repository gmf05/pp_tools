function data = make_seizure_pp(szX,patient_name,seizure_name,data_type);
  
  global DATA_DIR
  global MIN_REFRACT
  global SPIKE_THRESH_ECOG
  global SPIKE_THRESH_LFP
  global SPIKE_THRESH_MUA
  global DO_SAVE
  switch data_type
    case 'ECoG'
      thresh = SPIKE_THRESH_ECOG;
    case 'LFP'
      thresh = SPIKE_THRESH_LFP;
    case 'MUA'
      thresh = SPIKE_THRESH_MUA;
  end
  
  Name = [patient_name '_' seizure_name '_' data_type '_pp_thresh' num2str(thresh)];  
  d = szX.Data; % matrix of voltages (time = row, channel = col)
  N_channels = size(d,2);
  Fs = round(szX.SamplingRate);
  dt = 1/Fs;
  t = szX.Time; % NOTE: assumes time field exists
  
  % NOTE: assumes bad channels have already been filtered out with mod_sz
  % (1) preprocess with filters, etc.  
  d = preprocessing(d, t, data_type);

  % (2) find spikes channel by channel
  dn = 0*d';
  spikes = cell(1,N_channels);
  for i=1:N_channels
      spks = hilbertspike(d(:,i),thresh,MIN_REFRACT);      
      dn(i,spks) = 1;
      spikes{i} = t(spks);
  end
  
  % (3) remove any channels with very large/small spike counts
  temp = sum(dn'); % all spike counts
  cleantemp = removeoutliers(temp); % outliers removed
  out = setdiff(temp, cleantemp); % outliers
  N_out = length(out); 
  N_out = 0;
  out_ind = zeros(1,N_out);
%   for i=1:N_out, 
%     temp = find(temp==out(i));
%     Ni = length(temp);
%     out_ind(i:i+Ni-1) = temp;
%     i = i+Ni-1;
%   end
    
  good_ind = setdiff(1:N_channels,out_ind);
  dn = dn(good_ind,:);
  if isequal(class(szX.Labels),'char')
    szX.Labels = str2cell(szX.Labels);
  end
  Labels = {szX.Labels{good_ind}};
  fprintf(['\nRemoved ' num2str(N_out) ' ' data_type ...
    ' channels w/ outlying number of spikes.\n']);
  
  % (4) save, if desired
  data = pp_data(dn,t,Name,Labels);
  if DO_SAVE
    save([DATA_DIR '/' patient_name '/' patient_name ...
      '_' seizure_name '_' data_type '_spikes_thresh' num2str(thresh)],'spikes');
    save([DATA_DIR '/' patient_name '/' data.Name],'data','-v7.3');
  end
  
end

function d_post = preprocessing(d_pre, t, data_type)

dec = 1; % decimate, if necessary

if dec>1, d_post =[];
else d_post = d_pre; end;

%Get the sampling frequency.
Fs = 1/(t(2)-t(1));
fNQ = Fs/2;

switch data_type
    
    case 'ECoG'                        
        %Filters----
        %HighPass.
        f0 = [0   0.5   1   fNQ]/fNQ;
        z0 = [0   0   1     1]/fNQ;

        %LowPass.
        f1 = [1  119  119.5  fNQ]/fNQ;
        z1 = [1   1   0    0]/fNQ;

        %Bandstop
        f0S = [0  59.0   59.5  60.5    61 fNQ]/fNQ;
        z0S = [1     1     0    0       1   1]/fNQ;

        %Bandstop
        f1S = [0  119.0 119.5 120.5  121.0  fNQ]/fNQ;
        z1S = [1     1     0    0       1   1]/fNQ;

        b0L = firls(1000,f0,z0,'hilbert');
        b1L = firls(1000,f1,z1,'hilbert');
        b0S = firls(1500,f0S,z0S,'hilbert');
        b1S = firls(1500,f1S,z1S,'hilbert');

        for k=1:size(d_pre,2)
            if mod(k,10)==0, fprintf(['Channel #' num2str(k) '\n']), end        
            d0 = d_pre(:,k);
        %     d0 = decimate(d0, dec);           %Decimate (10).
            d0 = d0 - nanmean(d0);              %Subtract mean.
            d0 = filtfilt(b0L,1,d0);
            d0 = filtfilt(b1L,1,d0);
            d0 = filtfilt(b0S,1,d0);
            d0 = filtfilt(b1S,1,d0);
            d0 = zscore(d0);
            d_post(:,k) = d0;
        end
        
    case 'LFP'                        
        %Filters----
        %HighPass.
        f0 = [0   0.5   1   fNQ]/fNQ;
        z0 = [0   0   1     1]/fNQ;

        %LowPass.
        f1 = [1  199 199.5  fNQ]/fNQ;
        z1 = [1   1   0    0]/fNQ;

        %Bandstop
        f0S = [0  59.0   59.5  60.5    61 fNQ]/fNQ;
        z0S = [1     1     0    0       1   1]/fNQ;

        %Bandstop
        f1S = [0  119.0 119.5 120.5  121.0  fNQ]/fNQ;
        z1S = [1     1     0    0       1   1]/fNQ;

        %Bandstop
        f2S = [0  179.0 179.5 180.5  181.0  fNQ]/fNQ;
        z2S = [1     1     0    0       1   1]/fNQ;

        %Bandstop
        %f3S = [0  239.0 239.5 240.5  241.0  fNQ]/fNQ;
        %z3S = [1     1     0    0       1   1]/fNQ;

        b0L = firls(1000,f0,z0,'hilbert');
        b1L = firls(1000,f1,z1,'hilbert');
        b0S = firls(1500,f0S,z0S,'hilbert');
        b1S = firls(1500,f1S,z1S,'hilbert');
        b2S = firls(1500,f2S,z2S,'hilbert');
        %b3S = firls(1500,f3S,z3S,'hilbert'); % add'l notch filter

        for k=1:size(d_pre,2)
            if mod(k,10)==0, fprintf(['Channel #' num2str(k) '\n']), end        
            d0 = d_pre(:,k);
        %     d0 = decimate(d0, dec);           %Decimate (10).
            d0 = d0 - nanmean(d0);              %Subtract mean.
            d0 = filtfilt(b0L,1,d0);
            d0 = filtfilt(b1L,1,d0);
            d0 = filtfilt(b0S,1,d0);
            d0 = filtfilt(b1S,1,d0);
            d0 = filtfilt(b2S,1,d0);
            %d0 = filtfilt(b3S,1,d0);
            d0 = zscore(d0);            
            if dec>1, d0 = decimate(d0, dec); end;   %Decimate, if desired
            d_post(:,k) = d0;
        end
        % d_pre = d_pre0; clear d_pre0; % additional lines necessary if downsampling       
           
    case 'MUA'
        %Bandpass filters----
        % HighPass.
        fL = [0  299.5 300 fNQ]/fNQ;
        zL = [0   0   1     1]/fNQ;
        % LowPass.
        fH = [0  2999.5 3000 fNQ]/fNQ;
        zH = [0   0   1     1]/fNQ;

        bL = firls(1500,fL,zL,'hilbert');
        bH = firls(1500,fH,zH,'hilbert');

%         NT = 1:dec:size(d_pre,1); N_col = size(d_pre,2);
%         d_pre0 = zeros(NT,N_col);

        for k=1:size(d_pre,2)
            if mod(k,10)==0, fprintf(['Channel #' num2str(k) '\n']); end
            d0 = d_pre(:,k);
            d0 = d0 - mean(d0);
            d0 = filtfilt(bH,1,d0);
            d0 = filtfilt(bL,1,d0);
            d0 = zscore(d0);
            if dec>1, d0 = decimate(d0, dec); end;   %Decimate, if desired
            d_post(:,k) = d0;
        end
end

end

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
