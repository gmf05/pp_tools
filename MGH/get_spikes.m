function data = get_spikes(patient_name,seizure_name,data_type,thresh1,dT1,thresh2,dT2)
  
  if nargin<6, thresh2=thresh1; dT2=dT1; end
  
  global DATA; % path to data
  data_name = [patient_name '_' seizure_name '_' data_type '_' num2str(thresh1) '_' num2str(dT1) '_' num2str(thresh2)];
%   data_name = [patient_name '_' seizure_name '_' data_type '_' num2str(thresh1) '_' num2str(dT1) '_' num2str(thresh2) '_' num2str(dT2)];
  data_name0 = [patient_name ' ' seizure_name ' ' data_type ' @ thresh=' num2str(thresh1) ',' num2str(dT1) ',' num2str(thresh2) ',' num2str(dT2)];
  pp_filename = [DATA '/' patient_name '/' data_name '_pp.mat'];
  spikes_filename = [DATA '/' patient_name '/' data_name '_spikes.mat'];
  
  if exist(pp_filename,'file')
    fprintf(['Loading ' data_name0 '...']);
    load(pp_filename)
    fprintf('Done!\n');
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
      d = szX.Data'; % note: want channel x time 
%       d = zscore(d')'; % note: zscore expects time x channel
      
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
        spkind = spikefindA(-zscore(d(n,:)),thresh1,dT1,thresh2);
%         spkind = spikefindB(-zscore(d(n,:)),thresh1,dT1,thresh2,dT2);
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

function spikes = spikefindA(x,thresh1,dT,thresh2)
spikes=[];

ind = find(x>thresh1);
di = diff(ind);
jumps = [0 find(di>1)];
if length(jumps)==1, return, end
tOn = ind(jumps+1); tOn=tOn(1:end-1);
tOff = ind(jumps(2:end));

for j = 1:length(tOn);
  if tOff(j)-tOn(j)>dT
    [~,j0] = max(x(tOn(j):tOff(j)));
    j1 = tOn(j)+j0-1;
    if x(j1)>=thresh2
      spikes = [spikes tOn(j)+j0-1];
    end
  end
end

end


function spikes = spikefindB(x,thresh1,dT1,thresh2,dT2)
spikes=[];
T = length(x);
T0 = T-dT2;
t=dT1+1;

% time loop:
% march through time steps and check whether
% left-hand-slope > thresh1 & right-hand slope < thresh2 
% if so, local max is the spike
while t<T0  
  if x(t)-x(t-dT1)>=thresh1 & x(t)-x(t+dT2)>=thresh2
    % find local max:
    [~,t0] = max(x(t-dT1:t+dT2));
    t1 = t-dT1-1+t0;
    spikes = [spikes t1];
    t=t+dT2;
  else
    t=t+1; % no spike
  end
end

end