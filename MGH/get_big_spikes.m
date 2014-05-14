function data = get_big_spikes(patient_name,seizure_name,data_type)
  
  global DATA_DIR MIN_REFRACT SPIKE_THRESH_EEG SPIKE_THRESH_ECOG ...
    SPIKE_THRESH_LFP SPIKE_THRESH_MUA
  
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
  
  Name = [patient_name '_' seizure_name '_' data_type '_pp-big_thresh' num2str(thresh)];  
  pp_filename = [DATA_DIR '/' patient_name '/' Name '.mat'];
  
  if exist(pp_filename,'file')
    fprintf(['Loaded ' patient_name ' ' seizure_name ' @ thresh=' num2str(thresh) '\n']);
    load(pp_filename)
  else
    fprintf(['Cannot find spikes for ' patient_name ' ' seizure_name ' @ thresh = ' num2str(thresh) '\n']);
    fprintf('Loading raw data...\n');
    load([DATA_DIR '/' patient_name '/' patient_name '_' seizure_name '_LFP_ECoG_EEG.mat'],'sz');
    switch data_type
      case 'EEG'
        szX = sz.EEG;
        t = sz.ECoG.Time;
      case 'ECoG'
        szX = sz.ECoG;
        t = sz.ECoG.Time;
      case 'LFP'
        szX = sz.LFP;
        t = sz.LFP.Time;
      case 'MUA'
        sz.X = sz.LFP;
        t = sz.LFP.Time;
    end
    
    fprintf('Done!\nPreprocessing...\n');
    szX.Onset = sz.Onset; szX.Offset = sz.Offset;
    d = szX.Data;
    NT = size(d,1);
    N_channels = size(d,2);
    Fs = round(szX.SamplingRate);
    dt = 1/Fs;
    fNQ = Fs/2;
    
    spikes = cell(1,N_channels);
    start_ind = getclosest(t,szX.Onset);
    end_ind = getclosest(t,t(end)-szX.Offset);
    t = t(start_ind:end_ind) - t(start_ind);
    d = d(start_ind:end_ind,:);
    d = preprocessing(d, data_type);
    d = d'; % get in row = channel form
    save([DATA_DIR '/' patient_name '/' patient_name '_' seizure_name '_' data_type '_filtered'],'d','t');
    dn = 0*d;
    fprintf(['Done!\nFinding spikes...\n']);
    
    for n = 1:N_channels
      if mod(n,10)==1, fprintf(['Channel #' num2str(n) '\n']); end
      [spkind,amp] = hilbertspike(d(n,:),thresh,MIN_REFRACT);
      spikes{n} = t(spkind);
      amps{n} = amp;
      dn(n,spkind) = 1;
    end
    
    fprintf(['Done!\nSaving spikes...']);
    save([Name '_spikes.mat'],'spikes','amps','-v7.3');
    fprintf('Done!\n');

    % remove any channels with very large/small spike counts
    cumspks = sum(dn'); % all spike counts
    cleantemp = removeoutliers(cumspks); % outliers removed
    out = setdiff(cumspks, cleantemp); % set of outliers
    N_out = length(out);

    out_ind = [];
    count = 1;
    for i = 1:N_out
      ind_i = find(cumspks==out(i));
      Ni = length(ind_i);
      out_ind(count:count+Ni-1) = ind_i;
      count = count+Ni;
    end
    good_ind = setdiff(1:N_channels,out_ind);
    dn = dn(good_ind,:);
    if isequal(class(szX.Labels),'char')
      Labels = str2cell(szX.Labels(good_ind,:));
    elseif isequal(class(szX.Labels),'cell')
      Labels = {szX.Labels{good_ind}};
    end
    fprintf(['Removed ' num2str(N_out) ' ' data_type ...
      ' channels with too many/few spikes.\n']);
    
    data = pp_data(dn,t,Name,Labels);
    switch data_type
      case {'ECoG', 'EEG'}
        fprintf(['No downsampling.\n']);        
      case 'LFP', 
        fprintf(['Trying to downsample by 32x...']);
        try
          data = data.downsample(32); 
          fprintf('Done!\n');
        catch
          fprintf('Failed\n');
        end
      case 'MUA'
        fprintf(['Need to figure out whether to downsample here\n']);
    end
    
    fprintf(['Saving point process data object...']);
    save(pp_filename, '-v7.3','data');
    fprintf(['Done!\n\n']);
end

function d_post = preprocessing(d_pre, data_type)

switch data_type 

  %-------------------------- local field potential (LFP)
  case 'LFP'
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
  fH = [0   299.5   300.5   fNQ]/fNQ; zH = [0   0   1     1]/fNQ; % Highpass
  fL = [0 2999.5  3000.5  fNQ]/fNQ; zL = [1   1   0    0]/fNQ; % Lowpass
  bL = firls(2000,fH,zH);
  bH = firls(2000,fL,zL);
  d_post = d_pre;
  d_post = filtfilt(bH,1,d_post);
  d_post = filtfilt(bL,1,d_post);
  d_post = zscore(d_post);
  %--------------------------  multiunit activity (MUA)

  case 'ECoG'
  %-------------------------- filtered ECoG
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

  d_post = d_pre;
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

  d_post = d_pre;
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

