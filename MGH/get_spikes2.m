function data = get_spikes2(patient_name,seizure_name,data_type,thresh)
  
  global DATA;
 
  data_name = [patient_name '_' seizure_name '_' data_type '_thresh' num2str(thresh)];
  data_name0 = [patient_name ' ' seizure_name ' ' data_type ' @ thresh=' num2str(thresh)];
  pp_filename = [DATA '/' patient_name '/' data_name '_dpp.mat'];
  spikes_filename = [DATA '/' patient_name '/' data_name '_dspikes.mat'];
  filtered_filename = [DATA '/' patient_name '/' patient_name '_' seizure_name '_' data_type '_filtered_lo.mat'];
  
  if exist(pp_filename,'file')
    fprintf(['Loaded ' data_name0 '\n']);
    load(pp_filename)
  else    
    if exist(spikes_filename, 'file')
      fprintf(['Cannot find point process object for ' data_name0 '\n']);
      fprintf('Loading spike times...');
      load(spikes_filename)
      fprintf('Done!\n');
    else
      fprintf(['Cannot find spikes for ' data_name0 '\n']);
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
          case {'LFP','MUA'}
            szX = sz.LFP;
            t = sz.LFP.Time;
        end

        fprintf('Preprocessing...');
        szX.Onset = sz.Onset; szX.Offset = sz.Offset;
        labels = str2cell(szX.Labels);
        Fs = round(szX.SamplingRate);
        fNQ = Fs/2;
        d = szX.Data;
        
        % printing progress
        for i = 1:size(d,2), i, d(:,i) = preprocessing(d(:,i), data_type); end
        d = d'; % want channel x time
        save(filtered_filename,'-v7.3','d','t','labels');
        fprintf('Done!\n');
      end

      % get and save spikes, create raster      
      fprintf('Finding spikes...\n');
      N_channels = size(d,1);
      T = size(d,2);
      dn = 0*d;
      spikes = cell(1,N_channels);
      marks = cell(1,N_channels);
      if ~exist('Fs','var'), Fs = round(1/(t(2)-t(1))); end
      dW = round(.05*Fs); % how much shift to allow (# bins + or -)
      
      % loop over channels, finding spike indices for each
      for n = 1:N_channels
        if mod(n,10)==1, fprintf(['Channel #' num2str(n) '\n']); end
        spkind = derivspike(d(n,:),thresh);
        % post-process spkind to correspond to MINIMA over certain window
        % a spike s gets placed on the interval [s-dW, s+dW]          
        for i = 1:length(spkind)
          i0 = spkind(i) - dW;
          istart = max(1, spkind(i)-dW);
          iend = min(spkind(i)+dW, T);
          [~,shft] = min(d(n,istart:iend));
          spkind(i) = istart + shft - 1;
        end
        spkind = unique(spkind);
        % drop values 1 and 2 from spkind, if there
        spkind(spkind==1)=[]; spkind(spkind==2)=[];
        
        % save results
        dn(n,spkind) = 1;
        spikes{n} = t(spkind);
        d1 = diff(d(n,:)); d2 = diff(d1);        
        marks{n} = [d(n,spkind); d1(spkind-1); d2(spkind-2)];
      end

      fprintf('Done!\nSaving spikes...');
      save(spikes_filename,'-v7.3','spikes','marks','labels','t');
      fprintf('Done!\n');
    end
    
%     N_channels = length(spikes);
%     dn = zeros(N_channels,length(t));
%     for n = 1:N_channels      
%       dn(n,:) = hist(spikes{n},t);
%     end
    
    % save point process object
%     data = pp_data(dn,t);
    data = pp_data(dn,t,'name',data_name,'labels',labels,'marks',marks);
    data.labels = labels;
    data.name = data_name;
    data.marks = marks;
    
% %     % downsampling
% %     switch data_type
% %       case {'ECoG', 'EEG'}
% %         fprintf('No downsampling.\n');
% %       case 'LFP', 
% %         fprintf('Trying to downsample by 32x\n');
% %         try data = data.downsample(32); end;
% %       case 'MUA'
% %         fprintf('Need to figure out whether to downsample here\n');
% %     end
    
    fprintf('Saving point process data object...');
    save(pp_filename, '-v7.3','data');
    fprintf('Done!\n\n');
  end
end

function d_post = preprocessing(d_pre, data_type)

switch data_type 

  %-------------------------- local field potential (LFP)
  case 'LFP'
  fL = [0  59.5  60  fNQ]/fNQ; zL = [1   1   0    0]/fNQ; % Lowpass
  bL = firls(2000,fL,zL);
  d_post = d_pre;
  d_post = filtfilt(bL,1,d_post);
  %-------------------------- local field potential (LFP)
  
  %--------------------------  multiunit activity (MUA)
  case 'MUA'
  fH = [0  0.5     1   fNQ]/fNQ; zH = [0   0   1     1]/fNQ; % Highpass
  fL = [0 2999.5  3000.5  fNQ]/fNQ; zL = [1   1   0    0]/fNQ; % Lowpass
  bL = firls(2000,fH,zH);
  bH = firls(2000,fL,zL);
  d_post = d_pre;
  d_post = filtfilt(bH,1,d_post);
  d_post = filtfilt(bL,1,d_post);
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
  %-------------------------- filtered EEG
end
%   d_post = zscore(d_post);
end
end