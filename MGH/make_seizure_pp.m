function data = get_discharges(patient_name,seizure_name,data_type)
  
  global DATA_DIR
  global MIN_REFRACT
  global SPIKE_THRESH_EEG
  global SPIKE_THRESH_ECOG
  global SPIKE_THRESH_LFP
  global SPIKE_THRESH_MUA
  
  switch data_type
    case 'EEG'
      thresh = SPIKE_THRESH_EEG;
    case 'ECoG'
      thresh = SPIKE_THRESH_ECOG;
    case 'LFP'
      thresh = SPIKE_THRESH_LFP;
    case 'MUA'
      thresh = SPIKE_THRESH_MUA;
  end
  
  Name = [patient_name '_' seizure_name '_' data_type '_pp_thresh' num2str(thresh)];  
  pp_filename = [DATA_DIR '/' patient_name '/' Name '.mat'];
  
  d = szX.Data; % matrix of voltages (time = row, channel = col)  
  NT = size(d,1); 
  N_channels = size(d,2);
  Fs = round(szX.SamplingRate);
  dt = 1/Fs;
%   t = (1:NT)*dt; % time axis
  t = szX.Time;
  
  if exist(pp_filename,'file')
    load(pp_filename)
  else
    spikes = cell(1,N_channels);
    start_ind = getclosest(t,szX.Onset);
    end_ind = getclosest(t,t(end)-szX.Offset);
    t = t(start_ind:end_ind) - t(start_ind);
    d = d(start_ind:end_ind,:);
    d = preprocessing(d, t, data_type);
    dn = 0*d;
    for n = 1:N_channels, 
      spkind = hilbertspike(d,thresh,MIN_REFRACT);
      spikes{n} = t(spkind);
      dn(n,spkind) = 1;
    end 
    dn = dn(:,start_ind:end_ind);
    save([Name '_spikes.mat'],'spikes','-v7.3');

    % remove any channels with very large/small spike counts
    cumspks = sum(dn'); % all spike counts
    cleantemp = removeoutliers(cumspks); % outliers removed
    out = setdiff(cumspks, cleantemp); % outliers
    N_out = length(out); 
    out_ind = zeros(1,N_out);
    for i=1:N_out, 
      cumspks = find(cumspks==out(i));
      Ni = length(cumspks);
      out_ind(i:i+Ni-1) = cumspks;
      i = i+Ni-1;
    end
    good_ind = setdiff(1:N_channels,out_ind);
    dn = dn(good_ind,:);
    Labels = {szX.Labels{good_ind}};
    fprintf(['\nRemoved ' num2str(N_out) ' ' data_type ...
      ' channels w/ outlying number of spikes.\n']);
    
    data = pp_data(dn,t);
    switch data_type
      case {'ECoG', 'EEG'}
        fprintf(['No downsampling.\n']);        
      case 'LFP', 
        fprintf(['Downsampling by 32x\n']);
        data = data.downsample(32);
      case 'MUA'
        fprintf(['Need to figure out whether to downsample here\n']);
    end
    save(filename, '-v7.3','data');
  
end

function d_post = preprocessing(d_pre, data_type)

switch data_type 

  %-------------------------- local field potential (LFP)
  case 'LFP'
  Fs = 3e4; fNQ = Fs/2;
  fH = [0   299.5   300.5   fNQ]/fNQ; zH = [0   0   1     1]/fNQ; % Highpass
  fL = [1  200  200.5  fNQ]/fNQ; zL = [1   1   0    0]/fNQ; % Lowpass
  fS1 = [0  59.0   59.5  60.5    61 fNQ]/fNQ; zS1 = [1 1  0  0  1  1]/fNQ;  % Stop1
  fS2 = [0  119.0 119.5 120.5  121.0  fNQ]/fNQ; zS2 = [1  1  0  0  1 1]/fNQ; % Stop2
  fS3 = [0  179.0 179.5 180.5  181.0  fNQ]/fNQ; zS3 = [1  1  0  0  1 1]/fNQ; % Stop3
  bL = firls(2000,fH,zH);
  bH = firls(2000,fL,zL);
  bS1 = firls(2000,fS1,zS1);
  bS2 = firls(2000,fS2,zS2);
  bS3 = firls(2000,fS3,zS3);

  d_post = d_pre;
  d_post = filtfilt(bH,1,d_post);
  d_post = filtfilt(bL,1,d_post);
  d_post = filtfilt(bS1,1,d_post);
  d_post = filtfilt(bS2,1,d_post);
  d_post = filtfilt(bS3,1,d_post);
  d_post = zscore(d_post);
  %-------------------------- local field potential (LFP)

  case 'MUA'
  %--------------------------  multiunit activity (MUA)
  Fs = 3e4; fNQ = Fs/2;
  fH = [0   299.5   300.5   fNQ]/fNQ; zH = [0   0   1     1]/fNQ; % Highpass
  fL = [0 2999.5  3000.5  fNQ]/fNQ; zL = [1   1   0    0]/fNQ; % Lowpass
  bL = firls(2000,fH,zH);
  bH = firls(2000,fL,zL);
  d_post = v0;
  d_post = filtfilt(bH,1,d_post);
  d_post = filtfilt(bL,1,d_post);
  d_post = zscore(d_post);
  %--------------------------  multiunit activity (MUA)

  case 'ECoG'
  %-------------------------- filtered ECoG
  Fs = 5e3; fNQ = Fs/2;
  fH = [0  0.5   1.5   fNQ]/fNQ; zH = [0   0   1     1]/fNQ; % Highpass
  fL = [1  119.5 120.5  fNQ]/fNQ; zL = [1   1   0    0]/fNQ; % Lowpass
  fS1 = [0  59.0   59.5  60.5    61 fNQ]/fNQ; zS1 = [1 1  0  0  1  1]/fNQ;  % Stop1
  % fS2 = [0  119.0 119.5 120.5  121.0  fNQ]/fNQ; zS2 = [1  1  0  0  1 1]/fNQ; % Stop2
  % fS3 = [0  179.0 179.5 180.5  181.0  fNQ]/fNQ; zS3 = [1  1  0  0  1 1]/fNQ; % Stop3
  bL = firls(2000,fH,zH);
  bH = firls(2000,fL,zL);
  bS1 = firls(2000,fS1,zS1);
  % bS2 = firls(2000,fS2,zS2);
  % bS3 = firls(2000,fS3,zS3);

  d_post = v0;
  % d_post = d_post - nanmean(d_post);
  d_post = filtfilt(bH,1,d_post);
  d_post = filtfilt(bL,1,d_post);
  d_post = filtfilt(bS1,1,d_post);
  % d_post = filtfilt(bS2,1,d_post);
  % d_post = filtfilt(bS3,1,d_post);s
  d_post = zscore(d_post);
  %-------------------------- filtered ECoG

  case 'EEG'
  %-------------------------- filtered EEG
  Fs = 5e3; fNQ = Fs/2;
  fH = [0  0.5   1.5   fNQ]/fNQ; zH = [0   0   1     1]/fNQ; % Highpass
  fL = [1  119.5 120.5  fNQ]/fNQ; zL = [1   1   0    0]/fNQ; % Lowpass
  fS1 = [0  59.0   59.5  60.5    61 fNQ]/fNQ; zS1 = [1 1  0  0  1  1]/fNQ;  % Stop1
  % fS2 = [0  119.0 119.5 120.5  121.0  fNQ]/fNQ; zS2 = [1  1  0  0  1 1]/fNQ; % Stop2
  % fS3 = [0  179.0 179.5 180.5  181.0  fNQ]/fNQ; zS3 = [1  1  0  0  1 1]/fNQ; % Stop3
  bL = firls(2000,fH,zH);
  bH = firls(2000,fL,zL);
  bS1 = firls(2000,fS1,zS1);
  % bS2 = firls(2000,fS2,zS2);
  % bS3 = firls(2000,fS3,zS3);

  d_post = v0;
  % d_post = d_post - nanmean(d_post);
  d_post = filtfilt(bH,1,d_post);
  d_post = filtfilt(bL,1,d_post);
  d_post = filtfilt(bS1,1,d_post);
  % d_post = filtfilt(bS2,1,d_post);
  % d_post = filtfilt(bS3,1,d_post);s
  d_post = zscore(d_post);
  %-------------------------- filtered EEG
end
end

end

