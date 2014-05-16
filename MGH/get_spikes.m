function data = get_spikes(patient_name,seizure_name,data_type,thresh)
  
  global DATA; min_refract = 1;
 
  data_name = [patient_name '_' seizure_name '_' data_type '_thresh' num2str(thresh) '_pp'];
  pp_filename = [DATA '/' patient_name '/' data_name '.mat'];
  spikes_filename = [DATA '/' patient_name '/' patient_name '_' seizure_name '_' data_name '_thresh' num2str(thresh) '_spikes.mat'];
  filtered_filename = [DATA '/' patient_name '/' patient_name '_' seizure_name '_' data_type '_filtered.mat'];
  
  data_name0 = [patient_name ' ' seizure_name ' ' data_type ' @ thresh=' num2str(thresh)];
  
%   if exist(pp_filename,'file')
%     fprintf(['Loaded ' data_name0 '\n']);
%     load(pp_filename)
%   else    
    if exist(spikes_filename, 'file')
%       fprintf(['Cannot find point process object for ' data_name0 '\n']);
      fprintf('Loading spike times...');
      load(spikes_filename)
      fprintf('Done!\n');
    else
      fprintf(['Cannot find spikes for ' spike_request '\n']);
      if exist(filtered_filename,'file')
        fprintf('Loading processed data...');
        load(filtered_filename);
        fprintf('Done!\n');
      else
        fprintf('Cannot find processed data.\n');
        fprintf('Loading raw data...\n');
        load([DATA '/' patient_name '/' patient_name '_' seizure_name '_LFP_ECoG_EEG.mat'],'sz');
        fprintf('Done!\n');
        
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

        fprintf('Preprocessing...');
        szX.Onset = sz.Onset; szX.Offset = sz.Offset;
        labels = szX.Labels;
        Fs = round(szX.SamplingRate);
        fNQ = Fs/2;
        d = szX.Data;
        start_ind = getclosest(t,szX.Onset);
        end_ind = getclosest(t,t(end)-szX.Offset);
        t = t(start_ind:end_ind) - t(start_ind);
        d = d(start_ind:end_ind,:);
        
        % printing progress
        for i = 1:size(d,2), i, d(:,i) = preprocessing(d(:,i), data_type); end
%         d = preprocessing(d, data_type);
        d = d'; % want channel x time
        save(filtered_filename,'-v7.3','d','t','labels');
        fprintf('Done!\n');
      end

      % get and save spikes, create raster      
      fprintf('Finding spikes...\n');
      N_channels = size(d,1);
      spikes = cell(1,N_channels);
      amps = cell(1,N_channels);
      for n = 1:N_channels
        if mod(n,10)==1, fprintf(['Channel #' num2str(n) '\n']); end
        [spkind, amp] = hilbertspike(d(n,:),thresh,min_refract);
        spikes{n} = t(spkind);
        amps{n} = amp;
      end

      fprintf('Done!\nSaving spikes...');
      save(spikes_filename,'-v7.3','spikes','amps','labels','min_refract','t');
      fprintf('Done!\n');
    end

    % remove any channels with very large/small spike counts
    dn = zeros(length(spikes),length(t));
    for n = 1:length(spikes), dn(n,:) = hist(spikes{n},t); end
    cumspks = sum(dn,2); % all spike counts
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
    if isequal(class(labels),'char')
      labels = str2cell(labels(good_ind,:));
    elseif isequal(class(labels),'cell')
      labels = {labels{good_ind}};
    end
    fprintf(['Removed ' num2str(N_out) ' ' data_type ...
      ' channels with too many/few spikes.\n']);
    
    % save point process object
%     data = pp_data(dn,t);
    data = pp_data(dn,t,'name',data_name,'labels',labels,'marks',amps);
    data.labels = labels;
    data.name = data_name;
    data.marks = amps;
    
    switch data_type
      case {'ECoG', 'EEG'}
        fprintf('No downsampling.\n');
      case 'LFP', 
        fprintf('Trying to downsample by 32x\n');
        try data = data.downsample(32); end;
      case 'MUA'
        fprintf('Need to figure out whether to downsample here\n');
    end
    
%     fprintf('Saving point process data object...');
%     save(pp_filename, '-v7.3','data','data2');
%     fprintf('Done!\n\n');
%   end

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
  
  %--------------------------  multiunit activity (MUA)
  case 'MUA'
  fH = [0   299.5   300.5   fNQ]/fNQ; zH = [0   0   1     1]/fNQ; % Highpass
  fL = [0 2999.5  3000.5  fNQ]/fNQ; zL = [1   1   0    0]/fNQ; % Lowpass
  bL = firls(2000,fH,zH);
  bH = firls(2000,fL,zL);
  d_post = d_pre;
  d_post = filtfilt(bH,1,d_post);
  d_post = filtfilt(bL,1,d_post);
  d_post = zscore(d_post);
  %--------------------------  multiunit activity (MUA)
  
  %-------------------------- filtered ECoG
  case 'ECoG'
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

  %-------------------------- filtered EEG
  case 'EEG'
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