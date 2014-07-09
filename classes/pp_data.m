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
    T % number of time points
  end

  methods
    % Constructor
    function obj = pp_data(dn,t,varargin)
      
      obj.dn = dn;
      obj.N_channels = size(dn,1);
      obj.T = size(dn,2);
      if nargin<2, obj.t = 1:T;
      else obj.t = t; end
      obj.dt = obj.t(2) - obj.t(1);
      
      %
      % TO DO: parse varargin for auxillary data: marks, labels, etc
      %
      for n = 1:2:length(varargin)
        switch varargin{n}
          case 'name', obj.name = varargin{n+1};
          case 'labels', obj.labels = varargin{n+1};
          case 'marks',obj.marks = varargin{n+1};
        end
      end            
    end
    
    function obj = refresh(obj)
      obj = pp_data(obj.dn,obj.t,obj.name,obj.Labels);
    end
    
    % Returns a data object with the specified channels
    function obj2 = sub_data(obj,ind)      
      % check list ind to make sure entries are valid??
      dn = obj.dn(ind,:);
      if ~isempty(obj.labels), labels = {obj.labels{ind}}; else labels = {}; end
      obj2 = pp_data(dn,obj.t,'name',obj.name,'labels',labels);
      if ~isempty(obj.marks)
        obj2.marks = {obj.marks{ind}}; % for cell
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
      if obj.dt~=obj2.dt
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
    
    function plot(obj, plot_type, params)
      global PLOT_COLOR
      global FONT_SIZE 
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
          xlabel('time [s]')
          ylabel('rate [Hz]')
          title({ttl; plot_type});
        case 'psth'
          plot(obj.t,sum(obj.dn,1),PLOT_COLOR);
          xlabel('time [s]')
          ylabel('rate [Hz]')
          title({ttl ; ' - % channels active per dt'}); grid on;
        case 'heat'
          win_rates = zeros(obj.N_channels, N_windows);
          for i = 1:obj.N_channels
            count = 1;
            for n = 1:N_windows
              win_rates(i,n) = sum(obj.dn(i,count:count+win_bins-1));
              count = count+dW_bins;
            end
          end
          imagesc(win_rates);
          xlabel('time [s]')
          ylabel('channel')
          colorbar;
          title({ttl; 'firing rates'});
        case 'isi'
          if obj.dt>5e-4
            isi_axis = 1:10:ceil(1/obj.dt);
          else
            isi_axis = 1:100:ceil(.5/obj.dt);
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
          xlabel('time [s]');
          ylabel('ISI [ms]');
          title('inter-spike-interval (ISI) histogram');
          colorbar();
        case 'raster'          
          for i = 1:obj.N_channels
            ind = find(obj.dn(i,:));
            for k = 1:length(ind)
              plot([obj.t(ind(k)) obj.t(ind(k))], [i-0.5 i+0.5], PLOT_COLOR); hold on;
            end
          end
          xlabel('time [s]','fontsize',FONT_SIZE);
          ylabel('channel','fontsize',FONT_SIZE);
          title([ttl ' raster plot']);
          
        case 'raster-marks'
          min_mark = min([obj.marks{:}]);
          max_mark = max([obj.marks{:}]);
          Nsteps = 64;
          cstart = [1 0 0]; % red
          cend = [0 0 1]; % blue
          colors = interpolateColor(cstart, cend, Nsteps);
          for i = 1:obj.N_channels
            ind = find(obj.dn(i,:));
            for k = 1:length(ind)
              mk = obj.marks{i}(k);
              mkind = ceil((mk - min_mark)/(max_mark-min_mark)*Nsteps);
              col = colors(mkind,:);
              plot([obj.t(ind(k)) obj.t(ind(k))], [i-0.5 i+0.5], 'color', col); hold on;
            end
          end
          xlabel('time [s]','fontsize',FONT_SIZE);
          ylabel('channel','fontsize',FONT_SIZE);
          title([ttl ' raster plot']);
          
      end
      update_fig();
    end
    
    function obj0 = remove_outlier_counts(obj)
      counts = sum(obj.dn,2);
      [~,goodChan] = removeoutliers(counts);
      obj0 = obj.sub_data(goodChan);
    end

    function obj = downsample(obj, dT)
      fprintf(['Downsampling point process data by a factor of ' num2str(dT) '...']);
      t_old = obj.t;
      dn_old = obj.dn;
      
      obj.t = obj.t(1:dT:end);
      obj.dt = obj.dt*dT;
      obj.T = length(obj.t);      
      obj.dn = zeros(obj.N_channels,obj.T);
      
      for i = 1:obj.N_channels
        obj.dn(i,:) = hist(t_old(dn_old(i,:)>0), obj.t);
      end
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
      
    function spike_trigger_plot(obj,thresh,lockout)
      % 
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
      
% %       % against-all-spikes
%       figure(1);
%       global PLOT_COLOR
%       PLOT_COLOR = 'b'; obj.plot('raster'); PLOT_COLOR = 'r'; pause;

      % shifting trigger start and end by dL, dR 
      dL = 0.2; dR = 0.2; % [sec]
      dLbins = round(dL/obj.dt); dRbins = round(dR/obj.dt);
      
      istart=1;
      while istart < N_spks
        tstart = spkInfo(istart,1);
        iend=istart+1;
        while spkInfo(iend,1)-tstart<=lockout && iend < N_spks
          iend=iend+1;
        end        
        iend = iend-1;
        tend = spkInfo(iend,1);
        Nchan = length(unique(spkInfo(istart:iend,2)));
        if Nchan>=thresh
%           figure(2);
          obj.sub_time(tstart-dLbins:tend+dRbins).reset_time().plot('raster'); hold on;
%           figure(1);
%           obj.sub_time(tstart:tend).plot('raster'); hold on; % against-all-spikes
          pause();
        end
        istart = iend+1;
      end
    end
    
    function dcat = spike_trigger_plot2(obj,thresh,lockout)
      % 
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
      
% %       % against-all-spikes
      figure(1);
      global PLOT_COLOR
      PLOT_COLOR = 'b'; obj.plot('raster'); PLOT_COLOR = 'r'; pause;

      % shifting trigger start and end by dL, dR 
      dL = 0.2; dR = 0.2; % [sec]
      dLbins = round(dL/obj.dt); dRbins = round(dR/obj.dt);
      
      makeFlag = true;
      
      istart=1;
      while istart < N_spks
        tstart = spkInfo(istart,1);
        iend=istart+1;
        while spkInfo(iend,1)-tstart<=lockout && iend < N_spks
          iend=iend+1;
        end        
        iend = iend-1;
        tend = spkInfo(iend,1);
        Nchan = length(unique(spkInfo(istart:iend,2)));
        if Nchan>=thresh
%           figure(2);
          di = obj.sub_time(tstart-dLbins:tend+dRbins).reset_time();
          if makeFlag
            dcat = di;
            makeFlag = false;
          else
            dcat = dcat.concat(di);
          end
          
%           di.plot('raster'); hold on;
          figure(1);
          obj.sub_time(tstart:tend).plot('raster'); hold on; % against-all-spikes
          pause();
        end
        istart = iend+1;
      end
    end
    
    function obj0 = sort_mean_time(obj,thresh,lockout)
      % 
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
       
      istart=1;
      while istart < N_spks
        tstart = spkInfo(istart,1);
        iend=istart+1;
        while spkInfo(iend,1)-tstart<=lockout && iend < N_spks
          iend=iend+1;
        end
        iend = iend-1;
        tend = spkInfo(iend,1);
        Nchan = length(unique(spkInfo(istart:iend,2)));
        if Nchan>=thresh
          if ~exist('objcat','var')
            objcat = obj.sub_time(tstart:tend).reset_time();
          else
            objcat = objcat.concat(obj.sub_time(tstart:tend).reset_time());
          end
        end
        istart = iend+1;
      end
      
      srt = zeros(1,obj.N_channels);
      for cspk = 1:obj.N_channels
        srt(cspk) = mean(objcat.t(find(objcat.dn(cspk,:))));
      end
      [~,ord] = sort(srt);
      obj0 = obj.sub_data(ord);
    end
    
  end
end
