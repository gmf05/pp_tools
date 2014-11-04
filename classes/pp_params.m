classdef pp_params
    %   ---
    %   
    %   ---
    %   covariate_names
    %   covariate_knots
    %   covariate_bases:
    %   intrinsic_channel: channel being modeled
    %   ensemble_channels: lists of ensemble groups
    %   C: number of ensemble groups
    %   fit_method: 'glmfit', 'filt', or 'smooth'
    %   noise: noise parameters used for random walk in dynamic estimation
    %           (i.e. Kalman gain)
    %   window:
    %   downsample_est: factor for downsampling estimates (necessary for
    %                   filter, smoother)
    
  properties    
    covariate_names % e.g. 'rate', 'self-hist', 'ensemble'
    covariate_channels % channels corresponding to each covar
    covariate_knots % partition of time/lag axis
    covariate_bases % list of 'spline' or 'indicator' for each covar
    covariate_ind % lists of indices corresponding to each covar    
    response % channel being modeled
    fit_method % 'glmfit', 'filt', or 'smooth'
    link % link function for glmfit (default 'log')
    rs % how to rescale interspike intervals ('exp' default or 'identity')
    noise % noise (i.e. Kalman gain) parameters for filter/smoother
    s % tension parameter for spline interpolation (default s = 0.5)
    window % (default window = [1 0.2])
    downsample_est % factor for downsampling estimates (default 1)
    
  end
  
  methods
    function obj = pp_params(names,channels,knots,bases,fit_method,rs,noise)
      if nargin<1, names = {}; end
      if nargin<2, channels = {}; end
      if nargin<3, knots = {}; end
      if nargin<4, bases = {}; end      
      if nargin<5, fit_method = 'glmfit'; end      
      if isequal(fit_method, 'glmfit'), obj.link = 'log';
        else obj.link = ''; end
      if nargin<6, rs = 'exp'; end
      if nargin<7, noise=[]; end
      
      obj.covariate_names = names;
      obj.covariate_channels = channels;
      obj.covariate_knots = knots;
      obj.covariate_bases = bases;
      obj.covariate_ind = {};
      obj.fit_method = fit_method;
      obj.rs = rs;
      obj.noise = noise;
      obj.response = 1;
      
      N_covar_types = length(names);
      for n = 1:N_covar_types
        obj = obj.add_covar(names{n},channels{n},knots{n},bases{n});
      end
      
      obj.s = 0.5;
      obj.window = [1 0.2];
      obj.downsample_est = 1;
    end
    
    function obj = add_covar(obj, name, channels, knots, basis)
      if nargin<2, basis = 'indicator'; end;
      obj.covariate_names{end+1} = name;
      obj.covariate_channels{end+1} = channels;
      obj.covariate_knots{end+1} = knots;
      obj.covariate_bases{end+1} = basis;
      N = length(knots);
      N = N + 2*isequal(basis,'spline');
%       N = N - 1*isequal(basis,'indicator')*isequal(name,'rate');
      N = N - 1*isequal(basis,'indicator')*(channels==0);
      if isempty(obj.covariate_ind)        
        ind = 1:N;
      else
        ind = obj.covariate_ind{end}(end) + (1:N);
      end
      obj.covariate_ind{end+1} = ind;
    end
    
    function obj2 = get_covar(obj, i)
      obj2 = pp_params();
      obj2 = obj2.add_covar(obj.covariate_names{i},obj.covariate_channels{i},obj.covariate_knots{i},obj.covariate_bases{i});
      obj2.fit_method = obj.fit_method; 
    end
    
    function burn_in = get_burn_in(obj)
      burn_in = 0;
      for i = 1:length(obj.covariate_names)
        burn_in = max([burn_in, obj.covariate_knots{i}(end)]);
      end
    end

    function Xs = spline_Xi(obj, i)
      % i : covariate index -- makes a block of the spline matrix
      % like make_X_block
      knots = obj.covariate_knots{i};
      s = obj.s;
      
      if range(knots)>1, dtau = 1; else dtau = 0.01; end;
      tau = knots(1):dtau:knots(end);
      NT = length(tau);
      N = length(knots);
      onset = knots(1);
      offset = knots(end);
      intvls = diff(knots)*1/dtau;
      s_coeff = [-s  2-s s-2  s; 2*s s-3 3-2*s -s; ...
           -s   0   s   0;   0   1   0   0];
      
      if isequal(obj.covariate_bases{i}, 'indicator')
        if isequal(obj.covariate_names{i}, 'rate'), N = N-1; end
        Xs = eye(N); return;
      end
      
      Xs = zeros(NT,N+2);
      count=1;
      for n=1:N-1
        I = intvls(n); % length of interval (num. of bins)
        alphas = (0:I-1)./I;
        Xs(count:count+I-1, n+(0:3)) = [alphas'.^3 alphas'.^2 alphas' ones(I,1)] * s_coeff;
        count = count+I;
      end
      Xs(end, N-1:N+2) = [1 1 1 1] * s_coeff; % alpha = 1
    end


    function [tmins,tmaxs,N_windows] = window_taxis(obj,tmin,tmax)
        win = obj.window;
        my_case = length(win);
       
        switch my_case
          case 0
            tmins = tmin;
            tmaxs = tmax;
          case 1
            tmins = tmin : win : tmax - win;
            tmaxs = tmins + win;
          case 2
            win_size = win(1);
            dW = win(2);
            tmins = tmin:dW:tmax-win_size;
            tmaxs = tmins+win_size;            
        end
        N_windows = length(tmins);
    end   
    
    function [win_bins,dW_bins] = window_bins(obj,Fs)
      win_bins = round(obj.window(1)*Fs);
      dW_bins = round(obj.window(2)*Fs);
    end
  
  function Xs = splineX(obj,ind)
    s = obj.s;
    knots = obj.covariate_knots{ind};
    N = length(knots);
    s_coeff = [-s  2-s s-2  s; 2*s s-3 3-2*s -s; ...
           -s   0   s   0;   0   1   0   0];
    tau = knots(1):knots(end);
    NT = length(tau);
    intvls = diff(knots);
    Xs = zeros(NT,N+2);

    count=1;
    for n=1:N-1
      I = intvls(n); % length of interval (num. of bins)
      alphas = (0:I-1)./I;
      Xs(count:count+I-1, n+(0:3)) = [alphas'.^3 alphas'.^2 alphas' ones(I,1)] * s_coeff;
        count = count+I;
    end
    Xs(end, N-1:N+2) = [1 1 1 1] * s_coeff; % alpha = 1
  end
  
  function obj = set_covar_channels(obj,newchans)
    switch class(newchans)
      case 'double'
        for i = 1:length(newchans)
          obj.covariate_channels{i}=newchans(i);
        end
      case 'cell'
        for i = 1:length(newchans)
          obj.covariate_channels{i}=newchans{i};
        end
    end
  end
  
  end
  
end


