function data = get_spikes(patient_name,seizure_name,data_type,zthresh,dT)
  
  global DATA; % path to data 
  data_name = [patient_name '_' seizure_name '_' data_type '_' num2str(-zthresh) '_' num2str(dT)]; % note -zthresh
  data_name0 = [patient_name ' ' seizure_name ' ' data_type ' @ thresh=' num2str(zthresh) ',' num2str(dT)];
  pp_filename = [DATA '/' patient_name '/' data_name '_pp.mat'];
  spikes_filename = [DATA '/' patient_name '/' data_name '_spikes.mat'];
  
  if exist(pp_filename,'file')
    fprintf(['Loaded ' data_name0 '\n']);
    load(pp_filename)
  else
    if exist(spikes_filename, 'file')
      fprintf(['Cannot find point process object for ' data_name0 '\n']);
      fprintf('Loading spike times...');
      load(spikes_filename);      
      fprintf('Done!\n');
      
      N_channels = length(spikes);
      dn = zeros(N_channels,length(tfull));
      for n = 1:N_channels
        dn(n,:) = hist(spikes{n},tfull); 
      end
      
    else
      fprintf(['Cannot find spikes for ' data_name0 '\n']);
      fprintf('Loading raw data...');
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
      labels = str2cell(szX.Labels);
      d = szX.Data;
      d = zscore(d)'; % note: want channel x time 
      
      % get and save spikes, create raster      
      fprintf('Finding spikes...\n');
      dt = t(2)-t(1);     
      tfull = t(1):dt:t(end); % starting time axis t may have clipped intervals
      T = length(tfull); % here we swap "t" for the "full" time axis
      N_channels = size(d,1);
      dn = zeros(N_channels,T);
      spikes = cell(1,N_channels);
      marks = cell(1,N_channels);
      
      if ~exist('Fs','var'), Fs = round(1/(t(2)-t(1))); end
      
      % loop over channels, finding spike indices for each
      for n = 1:N_channels
        if mod(n,10)==1, fprintf(['Channel #' num2str(n) '\n']); end
        spkind = spikefind(zscore(d(n,:)),zthresh,dT);    
        spikes{n} = t(spkind);
        dn(n,:) = hist(spikes{n},tfull); % NOTE!!!: needed to modify "t" to be "full" time axis        
        marks{n} = [d(n,spkind)];
      end

      fprintf('Done!\nSaving spikes...');
      save(spikes_filename,'-v7.3','spikes','marks','labels','tfull');
      fprintf('Done!\n');
    end
       
    % save point process object 
    data = pp_data(dn,tfull,'name',data_name,'labels',labels,'marks',marks);
    data.labels = labels;
    data.name = data_name;
    data.marks = marks;
    
    fprintf('Saving point process data object...');
    save(pp_filename, '-v7.3','data');
    fprintf('Done!\n\n');
  end
end

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

t=dT+1;
while t<=T-dT
  t1 = t-dT;
  t2 = t+dT;
  dx1=x(t)-x(t1);
  dx2=x(t2)-x(t);
  if dx1*dx2<thresh && dx1<dx2    
    [~,ti] = min(x(t1:t2)); % find local min on [t1,t2]
    t0 = t1-1+ti; % index shifted from [t1,t2] to [1,T]
    spikes=[spikes t0];
    t=t2+dT;
  else
    t=t+1;
  end
end

end