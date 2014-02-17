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
  NT = size(d,1); 
  N_channels = size(d,2);
  Fs = round(szX.SamplingRate);
  dt = 1/Fs;
%   t = (1:NT)*dt; % time axis
  t = szX.Time;
  
  % NOTE: assumes bad channels have already been filtered out with mod_sz
  % (1) preprocess with filters, etc.  
  d = preprocessing(d, t, data_type);

  % (2) find spikes channel by channel
  dn = 0*d';
  for i=1:N_channels
      spks = hilbertspike(d(:,i),thresh,MIN_REFRACT);
      dn(i,spks) = 1;
  end
  
%   % (3) isolate seizure
%   start_ind = getclosest(t,szX.Onset);
%   end_ind = getclosest(t,t(end)-szX.Offset);
%   t = t(start_ind:end_ind) - t(start_ind); % normalize onset to t = 0
%   dn = dn(:,start_ind:end_ind);

  % (4) remove any channels with very large/small spike counts
  temp = sum(dn'); % all spike counts
  cleantemp = removeoutliers(temp); % outliers removed
  out = setdiff(temp, cleantemp); % outliers
  N_out = length(out); 
  out_ind = zeros(1,N_out);
  for i=1:N_out, 
    temp = find(temp==out(i));
    Ni = length(temp);
    out_ind(i:i+Ni-1) = temp;
    i = i+Ni-1;
  end
    
  good_ind = setdiff(1:N_channels,out_ind);
  dn = dn(good_ind,:);
  Labels = {szX.Labels{good_ind}};
  fprintf(['\nRemoved ' num2str(N_out) ' ' data_type ...
    ' channels w/ outlying number of spikes.\n']);
  
  % (5) save, if desired
  data = pp_data(dn,t,Name,Labels);
  if DO_SAVE,
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
        f1 = [1  200  205  fNQ]/fNQ;
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

