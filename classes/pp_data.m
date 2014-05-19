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
    marks % auxillary data for each spike, if desired
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
    function obj = sub_data(obj,ind)      
      % keep only valid indices:
      ind = intersect(1:obj.N_channels,ind);
      obj.dn = obj.dn(ind,:);
      obj.N_channels = size(obj.dn,1);   
      if ~isempty(obj.Labels)
        obj.Labels = {obj.Labels{ind}}; 
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

    function obj = concat(obj,obj2)
      % a couple general checks between old & new objects
      if obj.dt~=obj2.dt
        error('ERROR: data objects must have the same time resolution dt');
      elseif size(obj.dn,1)~=size(obj2.dn,1)
        error('ERROR: data objects must contain the same number of dimensions');
      end
      obj = pp_data([obj.dn obj2.dn],[obj.t obj.t(end) + (1:obj2.T)*obj2.dt]);
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
            plot(obj.t,obj.dn(i,:)+i-1,PLOT_COLOR); hold on;
          end
          xlabel('time [s]','fontsize',FONT_SIZE);
          ylabel('channel','fontsize',FONT_SIZE);
          title({ttl; 'raster plot'});
      end      
      update_fig();
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
    
  end
  
end
