classdef pp_data
% 
% pp_data: Point process data object
% Part of the Point Process Toolbox by Grant Fiddyment
% 
% Constructor: d = pp_data(dn, t, arguments)
% 
% dn: point process (i.e. event) data, being modeled (channels x time)
% t: time axis
% 
% --------- OPTIONAL ----------
% ARGUMENT    :   DESCRIPTION
% ----------------------------
% 'name'      :   data set title 
% 'labels'    :   identifier for each process
% 'marks'     :   continuous data associated with each process
% 

  properties
    name % title of the data set
    labels % identifier for each pp
    dn % point process data (rows = channels, cols = time)
    marks % auxillary data fors each spike, if desired
    N_channels % number of channels (rows)
    t % time axis
    dt % time resolution
    Fs % sampling frequency (sometimes easier than dt)
    T % number of time points
  end

  methods
    % Constructor
    function obj = pp_data(dn,t,varargin)
      
      obj.dn = dn;
      obj.N_channels = size(dn,1);      
      if nargin<2
        obj.T = size(dn,2);
        obj.t = 1:obj.T;
      else
        obj.t = t;
        obj.T = length(t);
      end
      
      if obj.T>1 % need more than one time point for sampling freq
        obj.dt = obj.t(2) - obj.t(1);
        obj.Fs = 1/obj.dt;
      end
     
      % parse varargin for name, marks, labels
      for n = 1:2:length(varargin)
        switch varargin{n}
          case 'name', obj.name = varargin{n+1};
          case 'labels', obj.labels = varargin{n+1};
          case 'marks',obj.marks = varargin{n+1};
        end
      end            
    end
    
    function obj = refresh(obj)
      obj = pp_data(obj.dn,obj.t,'name',obj.name,'labels',obj.labels,'marks',obj.marks);
    end
    
    % Returns a data object with the specified channels
    function obj2 = sub_data(obj,ind)
      % check list ind to make sure entries are valid??
      dn = obj.dn(ind,:);
      if ~isempty(obj.labels), labels = {obj.labels{ind}}; else labels = {}; end
      obj2 = pp_data(dn,obj.t,'name',obj.name,'labels',labels);
      if ~isempty(obj.marks)
        obj2.marks = reshape({obj.marks{:,ind}},[],obj2.N_channels); % for cell
%         obj2.marks = obj.marks(ind,:); % for matrix
      end
    end
    
    function obj = sub_time(obj,varargin)
      if length(varargin)==1
        ind = varargin{1};
      else
        beg_ind = getclosest(obj.t,varargin{1});
        end_ind = getclosest(obj.t,varargin{2});
        ind = beg_ind:end_ind;
      end
      
      % get marks within particular time interval
      for n = 1:obj.N_channels
        beg_mark = sum(obj.dn(n,1:ind(1)-1))+1;
        end_mark = sum(obj.dn(n,1:ind(end)));
        obj.marks{n} = obj.marks{n}(:,beg_mark:end_mark);
      end
      
      obj.dn = obj.dn(:,ind);
      obj.t = obj.t(ind);
      obj.T = length(obj.t);
      
    end
    
    function obj = sub_time_fast(obj,varargin)
      % keeps all marks rather than distinguishing
      % which are related to the new time window
      % this approach is meant for quick and dirty
      % applications
      obj.marks = {};
      if length(varargin)==1
        ind = varargin{1};
      else
        beg_ind = getclosest(obj.t,varargin{1});
        end_ind = getclosest(obj.t,varargin{2});
        ind = beg_ind:end_ind;
      end
      
      obj.dn = obj.dn(:,ind);
      obj.t = obj.t(ind);
      obj.T = length(obj.t);
    end

    function [t,T,dt,tmin,tmax] = get_time(obj)
      t=obj.t;
      T=obj.T;
      dt=obj.dt;
      tmin=t(1)-dt;
      tmax=t(end);
    end
    
    function obj = reset_time(obj)
      obj.t = (1:obj.T)*obj.dt;
    end

    function obj = concat(obj,obj2)
      % a couple general checks between old & new objects
      if abs(obj.dt-obj2.dt)>1e-10
        error('ERROR: data objects must have the same time resolution dt');
      elseif size(obj.dn,1)~=size(obj2.dn,1)
        error('ERROR: data objects must contain the same number of dimensions');
      end
      obj.dn = [obj.dn obj2.dn];
      obj.t = [obj.t obj2.t];
      obj.marks = [obj.marks obj2.marks];
      obj.T = length(obj.t);
      % check labels to ensure they match??
    end
    
    function obj = add_channels(obj,obj2)
      obj.dn = [obj.dn; obj2.dn];
      obj.N_channels = obj.N_channels+size(obj2.dn,1);
      try
        obj.labels = {obj.labels{:} obj2.labels{:}};
      catch
        0;
      end
    end
    
    function plot(obj, plot_type, params)
      global PLOT_COLOR
      if nargin<2, plot_type = 'rate'; end
      if nargin<3, params = pp_params(); end
      
      [tmins,~,N_windows] = params.window_taxis(obj.t(1),obj.t(end));
      %%% HACK: something is broken here
      N_windows = N_windows - 1;  tmins = tmins(1:N_windows);
      %%%
      win_bins = floor(params.window(1) / obj.dt);
      dW_bins = floor(params.window(2) / obj.dt);      
      if ~isempty(obj.name)
        ttl = regexprep(obj.name,'_',' ');
        ttl = regexprep(ttl, ' pp thresh', '');
      else
        ttl = '';
      end
      
      win_rate = zeros(obj.N_channels,N_windows);
      switch plot_type
        case 'rate'
          for i = 1:obj.N_channels
            count = 1;
            for n = 1:N_windows
              win_rate(i,n) = sum(obj.dn(i,count:count+win_bins-1));
              count = count+dW_bins;
            end
            plot(tmins,win_rate,PLOT_COLOR); hold on; grid on;
            
          end
          xlabel('Time [s]')
          ylabel('rate [Hz]')
          title({ttl; plot_type});
        case 'psth'
          plot(obj.t,sum(obj.dn,1),PLOT_COLOR);
          xlabel('Time [s]')
          ylabel('# Channels')
          title({ttl ; ' - no. channels active per dt'}); grid on;
        case 'heat'
          Ts = obj.t(1:dW_bins:end-win_bins);
          win_rates = zeros(obj.N_channels, N_windows);
          for i = 1:obj.N_channels
            count = 1;
            for n = 1:N_windows
              win_rates(i,n) = sum(obj.dn(i,count:count+win_bins-1));
              count = count+dW_bins;
            end
          end
          imagesc(Ts,1:obj.N_channels,win_rates);
          xlabel('Time [s]')
          ylabel('Channel')
          colorbar;
          title({ttl; 'firing rates'});
        case 'isi'
          isi_axis = 1:50:2000;
          isi_list=[];
          for n = 1:obj.N_channels, isi_list=[isi_list round(diff(find(obj.dn(n,:)))*obj.dt*1e3)]; end
          [isi_hist,isi_x] = hist(isi_list,isi_axis);
          bar(isi_x,isi_hist);
          xlabel('ISI [ms]');
          ylabel('Count');
        case 'isi-heat'
          if obj.dt>5e-4
            isi_axis = 1:10:ceil(1/obj.dt);
          else
            isi_axis = 1:100:ceil(1/obj.dt);
%             isi_axis = 1:100:ceil(2/obj.dt);
          end
          isi_hist = zeros(length(isi_axis),N_windows);
          count = 1;
          for n = 1:N_windows
            % get histogram of ISI over all channels
            for i = 1:obj.N_channels
              isi_i_n = diff(find(obj.dn(i,count:count+win_bins-1)));
              % drop first spike since ISI may be affected by window
              isi_i_n = isi_i_n(2:end);
              [A,B] = hist(isi_i_n,isi_axis);              
              isi_hist(:,n) = isi_hist(:,n) + A';
            end
            isi_hist(:,n) = isi_hist(:,n) ./ sum(isi_hist(:,n));
            count = count+dW_bins;
          end
          imagesc(tmins,round(isi_axis*obj.dt*1e3),isi_hist);
          xlabel('Time [s]');
          ylabel('ISI [ms]');
          title('inter-spike-interval (ISI) histogram');
          colorbar();
        case 'raster'
          gca(); hold on;
          % fixing plotting so it's not horrendously slow...
          for i = 1:obj.N_channels
            ind = find(obj.dn(i,:));
            plot(obj.t(ind), i*ones(1,length(ind)), [PLOT_COLOR '.']);
          end
          ylim([0.5, obj.N_channels+0.5]);
          xlim([obj.t(1) obj.t(end)]);
          xlabel('Time [s]');
          ylabel('Channel');
          title([ttl ' Raster Plot']);
          axis ij;
        case 'raster2'
          gca(); hold on;
          % old plotting routine: connected spike trains but SLOW
          % -- OK if only a few channels / small time interval
          for i = 1:obj.N_channels
            ind = find(obj.dn(i,:));
            for k = 1:length(ind)
              plot([obj.t(ind(k)) obj.t(ind(k))], [i-0.5 i+0.5], 'Color', PLOT_COLOR);
            end
          end
          ylim([0.5, obj.N_channels+0.5]);
          xlim([obj.t(1) obj.t(end)]);
          xlabel('Time [s]');
          ylabel('Channel');
          title([ttl ' Raster Plot']);
          axis ij;
          
          
          
        case 'raster-marks'
          j=3; % row of mark process
%           min_mark = min(min([obj.marks{:}]));
%           max_mark = max(max([obj.marks{:}]));
          Nsteps = 64;
%           cstart = [1 0 0]; % red
%           cend = [0 0 1]; % blue
%           colors = interpolateColor(cstart, cend, Nsteps);
%           colors = colormap(jet);
          colors = colormap(hsv);
          
          for i = 1:obj.N_channels
            ind = find(obj.dn(i,:));
            marks_norm = normalize(obj.marks{i}(j,:));
            col_ind = round(marks_norm .* (Nsteps-1)) + 1;
            for k = 1:length(ind)
              plot([obj.t(ind(k)) obj.t(ind(k))], [i-0.5 i+0.5], 'color', colors(col_ind(k),:),'linewidth',3.5);
              hold on;
            end
          end
          xlabel('Time [s]');
          ylabel('Channel');
          title([ttl ' Raster Plot']);
          axis ij
      end
    end
    
    function obj0 = remove_outlier_counts(obj)
      counts = sum(obj.dn,2);
      [~,goodChan] = removeoutliers(counts);
      obj0 = obj.sub_data(goodChan);
    end

    function obj = downsample(obj, dT)
      fprintf(['Downsampling point process data by a factor of ' num2str(dT) '...']);
%       t_old = obj.t;
%       dn_old = obj.dn;
      obj.t = obj.t(1:dT:end);
      obj.dt = obj.dt*dT;
      obj.Fs = 1/obj.dt;
      obj.T = length(obj.t);     
      
      % method 1
%       obj.dn = zeros(obj.N_channels,obj.T);      
%       for i = 1:obj.N_channels
%         obj.dn(i,:) = hist(t_old(dn_old(i,:)>0), obj.t);
%       end
           
%       % method 2
      obj.dn = cumdownsample(obj.dn,dT);
      
      fprintf('Done!\n');
      
      if any(obj.dn(:)>1)
        fprintf(' \nWARNING: Downsampling has caused values > 1 \n');
      end            
    end

    function obj = upsample(dT,method,delta)
      if nargin<3, method = 'standard'; end
      if nargin<4, delta = []; end
      obj.dt = obj.dt/dT;
      obj.t = obj.t(1):obj.dt:obj.t(end);
      obj.T = length(obj.t);
      N = size(obj.dn,1);
      
      dn_old = obj.dn;
      obj.dn = zeros(N,obj.T);
      
      for n=1:N
        spks = obj.t(dn_old(n,:)>0);
        % adds jitter, if desired
        switch method
          case 'jitter'
            for s = 1:length(spks)
                spks(s) = spks(s) + delta*rand;
            end
        end
        obj.dn(n,:) = hist(spks,obj.t);     
      end
    end

    function obj0 = jitter(obj,delta,N_samples,method,R)
      rng('shuffle');
      % sub_data.jitter(channel,delta,N_samples,method,R)
      % NOTE: this method assumes the invoking object obj is
      % a sub-data object, i.e. d.sub_data(#), where d is a larger data
      % object and # is the desired channel number. If the object d
      % has only 1 channel this is unnecessary.
      %
      % delta = size of jitter window
      % N_samples = number of draws to make when resampling
      % method = 'basic', 'pattern', or 'interval'
      % R = window where history-dependent structure is preserved (if
      %           method = 'pattern')
      %
      % 
      %
      if nargin<3, N_samples = 1; end
      if nargin<4, method = 'basic'; end
      if nargin<5 && isequal(method,'pattern'), error('Musxt specify pattern length R'); end;

      obj0 = obj;
      dn_jitter = zeros(N_samples,obj0.T);
      % here we follow the pattern jitter & interval jitter algorithms
      % of Harrison & Geman (2008):

      spk_ind = find(obj0.dn(1,:)>0);
      spk_ind_rs = zeros(length(spk_ind),N_samples); % resampled spike indices

      switch method
        case 'basic' 
        % jitter all spikes w/ uniform [0,delta] noise, don't preserve patterns
        % also, doesn't center jittered samples around originals
          for i = 1:length(spk_ind)
            spk_ind_rs(i,:) = spk_ind(i) + floor(rand(1,N_samples)*(delta+1));
          end
      case 'pattern'                
        spk_ind_rs(1,:) = spk_ind(1); % intialize resampled spikes
        spk_ind_rs(end,:) = spk_ind(end);
        for i = 2:length(spk_ind)
          dX = spk_ind(i)-spk_ind(i-1);
          if dX<R, X(i,:) = X(i-1,:)+dX;
          else
            % compute hi
            hi = 0;
            % compute Bi
            Bi = Bi_end;
            for j = length(spk_ind):-1:i
              Bi = 0; 
            end
            % compute pi
            pi = Bi ./ sum(Bi);
            % draw re-samples
            spk_ind_rs(i,:) = 0;
          end
        end
      case 'interval'
          [];
      end

      for n = 1:N_samples, dn_jitter(n,spk_ind_rs(:,n)') = 1; end            
      obj0.dn = dn_jitter;
    end
    
    
    function spkInfo = raster_ind(obj)
      psth = sum(obj.dn);
      spkInd = find(psth);
      N_spks = sum(psth);
      spkInfo = zeros(N_spks,2);
      count=0;
      for n = spkInd
        temp = find(obj.dn(:,n));        
        N = length(temp);
        spkInfo(count+(1:N),1) = n;
        spkInfo(count+(1:N),2) = temp;
        count = count+N;
      end
    end
    
    function intvls = spike_trigger(obj,thresh,lockout,dL,dR)
      % set parameters
      if nargin<4, dL=0; end
      if nargin<5, dR=0; end
      dLbins = round(dL/obj.dt); dRbins = round(dR/obj.dt);

      % get array of spike indices, channels
      spkInfo = obj.raster_ind();
      Nspks = size(spkInfo,1);

      % iterate over spikes, finding windows where enough
      % channels spike
      intvls = [];
      spk1=1;
      while spk1 < Nspks
        istart = spkInfo(spk1,1);
        % find last spike before lockout period & channel repeat
        spk2 = spk1+1;
        Nchans = length(unique(spkInfo(spk1:spk2,2)));
        while spkInfo(spk2,1)-istart<=lockout && spk2 < Nspks && Nchans==spk2-spk1+1
          spk2 = spk2+1;
          chans = spkInfo(spk1:spk2,2);
          uchans = unique(chans);
          Nchans = length(uchans);
        end
        spk2 = spk2-1;

        if Nchans>=thresh
          iend = spkInfo(spk2,1);
          intvls = [intvls; istart-dLbins min(iend+dRbins,obj.T)];
        end

        spk1 = spk2+1;
      end
      
      % make sure detected intervals lie within data's time axis
      if intvls(1) < 1, intvls(1) = 1; end
      if intvls(end) > obj.T, intvls(end) = obj.T; end
    end
    
    function intvls = spike_trigger2(obj,thresh,lockout,dL,dR)
      dLbins = round(dL/obj.dt); dRbins = round(dR/obj.dt);

      % get array of spike indices, channels
      spkInfo = obj.raster_ind();
      Nspks = size(spkInfo,1);

      % iterate over spikes, finding windows where enough
      % channels spike
      intvls = [];
      spk1=1;
      while spk1 < Nspks
        hasRepeat = false; % prevent channel from appearing twice in one wave
        istart = spkInfo(spk1,1);
        % find last spike before delay of 'lockout' occurs, repeat
        spk2 = spk1;
        while spk2 < Nspks && spkInfo(spk2+1,1) - spkInfo(spk2,1) <= lockout && ~hasRepeat
          spk2 = spk2+1;
          counts = histc(spkInfo(spk1:spk2,2),1:obj.N_channels);
          hasRepeat = sum(counts>1);
        end
        if hasRepeat, spk2 = spk2-1; end
        Nchans = length(unique(spkInfo(spk1:spk2,2)));

        % ensure at least 'thresh' unique channels spike
        if Nchans>=thresh
          iend = spkInfo(spk2,1);
          intvls = [intvls; istart-dLbins min(iend+dRbins,obj.T)]; %#ok
        end
        spk1 = spk2+1;
      end
      % make sure detected intervals lie within data's time axis
      if isempty(intvls), return; end
      intvls(intvls<1) = 1;
      intvls(intvls>obj.T) = obj.T;
    end
    
    function obj0 = spike_trigger_hist(obj,intvls,time0)
      objcat = obj.sub_time_fast(intvls(1,1):intvls(1,2)).reset_time();
      objcat.t = objcat.t - time0(1);
      for i = 2:size(intvls,1)
        obji = obj.sub_time_fast(intvls(i,1):intvls(i,2)).reset_time();
        obji.t = obji.t + time0(1);
        objcat = objcat.concat(obji);
      end
      
      obj0 = obj;
      obj0.dn = zeros(obj.N_channels, length(time0));
      for n = 1:obj.N_channels
        obj0.dn(n,:) = hist(objcat.t(objcat.dn(n,:)>0), time0);
      end
      
    end
    
    function mov = spike_trigger_plot(obj,intvls,dTL)
      N = size(intvls,1);
      global PLOT_COLOR
      rgb = jet(N);
      for n = 1:N
%         obj.sub_time_fast(intvls(i,1):intvls(i,2)).reset_time().plot('raster');
        obj0 = obj.sub_time_fast(intvls(n,1):intvls(n,2)).reset_time();
        obj0.t = obj0.t - dTL;
        PLOT_COLOR = rgb(n,:);
        obj0.plot('raster2')
        hold on;
%         pause(0.05);
%         mov(n) = getframe();
      end
    end
        
    function [obj0, ord] = sort_mean_time(obj,intvls)
%       intvls = obj.spike_trigger(thresh,lockout);
%       intvls(end,:)=[]; % weird bug, last interval is no good??
      objcat = obj.sub_time_fast(intvls(1,1):intvls(1,2)).reset_time();
      for i = 2:size(intvls,1)
        objcat = objcat.concat(obj.sub_time_fast(intvls(i,1):intvls(i,2)).reset_time());
      end
      
      srt = zeros(1,obj.N_channels);
      for cspk = 1:obj.N_channels
        srt(cspk) = mean(objcat.t(find(objcat.dn(cspk,:))));
      end
      [~,ord] = sort(srt);
      obj0 = obj.sub_data(ord);
    end
    
    function [means, stds] = trigger_stats(obj, intvls)
      spkInfo = obj.raster_ind();
      spkInfo0 = [spkInfo 0*spkInfo(:,1)];
      j1 = 1;
      % for each wave, which spikes belong to it?
      for k = 1:size(intvls,1)
        while spkInfo(j1,1) < intvls(k,1)
          j1 = j1+1;
        end
        j2 = j1;
        while spkInfo(j2,1) < intvls(k,2) 
          j2 = j2+1;
        end
        if spkInfo(j2,1) > intvls(k,2)
          j2 = j2-1;
        end
        spkInfo0(j1:j2, 3) = k;
        j1 = j2+1;
      end
      [~, ind] = sort(spkInfo(:,2));
      spkInfo0 = spkInfo0(ind, :);

      means = zeros(1, obj.N_channels);
      stds = zeros(1, obj.N_channels);
      for i = 1:obj.N_channels
        ind = find(spkInfo0(:, 2)==i); % get spikes from channel i
        ind0 = ind(spkInfo0(ind, 3)>0); % get spikes that fall in waves
        N = length(ind0); % number of spikes falling within waves
        starts = intvls(spkInfo0(ind0, 3), 1);
        delays = (spkInfo0(ind0, 1) - starts) * obj.dt * 1e3; % time since wave start [in ms]
        means(i) = nanmean(delays);
        stds(i) = nanstd(delays);
%         stds(i) = quantile(delays, 0.75) - quantile(delays, 0.25);
      end
%       error('test') % debug stop
    end
  end
end
