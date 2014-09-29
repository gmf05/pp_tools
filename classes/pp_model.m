  classdef pp_model
%
%   properties
%     b % parameters (column vector)
%     W % covariance matrix
%     X % design (data) matrix
%     y % response process
%     CIF % conditional intensity function 
%     fit_method % glmfit, filt, or smooth
%     link % link function name, if using glmfit 
%     stats % structure like the one returned by glmfit
%     LL % log-likelihood
%     dev % deviance 
%     AIC % Akaike information criterion (penalized LL)    
%     rsISI % rescaled inter-spike-intervals ("waiting times")
%     KS % Kolmogorov-Smirnov statistic (are rescaled ISI uniform?)
%   end
%
%   methods
%     m = m.fit(d,p)    
%     Estimates model parameters & covariance
%     INPUT: point process data (d), parameters (p)
%     OUPUT: estimated point process model (m)
%
%     m.gof_plot(d)
%     Shows goodness-of-fit summary for model m
%     INPUT: point process data (d)
%
%     m.plot(d,p)
%     Plots model estimates
%     INPUT: point process data (d), parameters (p)
%     
%   end

  properties
    b % parameters (column vector)
    W % covariance matrix
    X % design (data) matrix
    y % response process
    CIF % conditional intensity function 
    fit_method % glmfit, filt, or smooth
    link % link function name, if using glmfit 
    stats % structure like the one returned by glmfit
    LL % log-likelihood
    dev % deviance 
    AIC % Akaike information criterion (penalized LL)    
    rs % how to rescale inter-spike-intervals ('exp' default or 'identity')
    rsISI % rescaled inter-spike-intervals ("waiting times")
    KS % Kolmogorov-Smirnov statistic (are rescaled ISI uniform?)
  end
    
  methods
    % Constructor
    function obj = pp_model(b,W,X,y,C,link)      
      obj.rs = 'exp';
      if nargin>0
        obj.b = b;
        obj.W = W;
        if nargin>2
          obj.X = X;
          obj.y = y;
          obj.CIF = C;
          obj.link = link;
          obj = obj.calcGOF();
        end
      end
    end
        
    function obj = fit(obj, d, p)
    % pp_model.fit(d, p)
    % INPUTS:
    % d -- point process data object
    % p -- point process params object
    %
    %
      warning(''); % clear last warning 
      fprintf(['\nFitting point process model...\n']);      
      obj.fit_method = p.fit_method;
      obj.link = p.link;
      obj.rs = p.rs;
      N_cov_types = length(p.covariate_names);
      obj.y = d.dn(p.response,:)'; % response variable      
      
      % make design matrix
      N_cov = p.covariate_ind{end}(end); % total # of covariates
      fs_update_ind = p.covariate_ind{1}+1:N_cov;
      obj.X = ones(d.T,N_cov);
      fprintf(['Building design matrix...\n']);
      obj = obj.makeX(d,p);      
      fprintf(['Done!\n']);
      
      % trim burn-in period
      burn_in = p.get_burn_in();
      obj.X = obj.X(burn_in+1:end,:);
      obj.y = obj.y(burn_in+1:end);            
      
      fprintf(['Estimating parameters...\n']);
%       % MATLAB glmfit routine:
%       [b,dev,stats] = glmfit(obj.X,obj.y,'poisson','link',obj.link,'constant','off');

      % custom glmfit routine:
      [b,stats] = obj.glmfit0(obj.X,obj.y,obj.link);

      switch obj.fit_method
        case 'glmfit'
          obj.b = b;
          obj.W = stats.covb;
%           obj.CIF = glmval(b,obj.X,obj.link,'constant','off');
          obj.CIF = exp(obj.X*b);
          obj.stats=stats;
        
        case 'filt'                    
          % initialize arrays
          NT = length(burn_in:p.downsample_est:d.T);
          bs = cell(1,NT);
          Ws = cell(1,NT);
          obj.CIF = zeros(d.T - burn_in, 1);
          X_bwd = flipud(obj.X);
          y_bwd = fliplr(obj.y);
%           [b_bwd,~,stats_bwd] = glmfit(X_bwd,y_bwd,'poisson','link',obj.link,'constant','off');
          [b_bwd,stats_bwd] = obj.glmfit0(X_bwd,y_bwd,obj.link);
          W_bwd = stats_bwd.covb;
          
          % noise (Kalman gain) matrix:
          noise_mtx = zeros(N_cov);
          for n = 1:N_cov_types
            for i = p.covariate_ind{n}
              noise_mtx(i, i) = p.noise(n);
            end
          end
          noise_mtx = noise_mtx(fs_update_ind,fs_update_ind);
          
          fprintf(['Filtering backward...']);
          for t = 1:d.T-burn_in
            dL_dB = X_bwd(t,fs_update_ind);
            lt = exp(X_bwd(t,:)*b_bwd);
            A = W_bwd(fs_update_ind,fs_update_ind)+noise_mtx;
            U = dL_dB';
            V = dL_dB;
            C0 = lt;
            C1 = (1/C0+V*A*U);
            % matrix inversion lemma
            W_bwd(fs_update_ind,fs_update_ind) = A-A*U*(C1\V)*A;
            dB = W_bwd(fs_update_ind,fs_update_ind)*dL_dB'*(obj.y(t) - lt);
            b_bwd(fs_update_ind) = b_bwd(fs_update_ind) + dB;
          end
%           b = b_bwd;
%           b = rand(size(b));
          W = W_bwd;
          fprintf(['Done!\n']);
          
          % set initial values          
          bs{1} = b;
          Ws{1} = W;
          count = 2;
          
          fprintf(['Filtering forward...']);
          for t = 1:d.T-burn_in
            dL_dB = obj.X(t,fs_update_ind);
            lt = exp(obj.X(t,:)*b);
            A = W(fs_update_ind,fs_update_ind)+noise_mtx;
            U = dL_dB';
            V = dL_dB;
            C0 = lt;
            C1 = (1/C0+V*A*U);
            % matrix inversion lemma:
            W(fs_update_ind,fs_update_ind) = A-A*U*(C1\V)*A;
            dB = W(fs_update_ind,fs_update_ind)*dL_dB'*(obj.y(t) - lt);
            b(fs_update_ind) = b(fs_update_ind) + dB;
            obj.CIF(t) = exp(obj.X(t,:)*b);
            
            if mod(t, p.downsample_est)==0
              bs{count} = b;
              Ws{count} = W;
              count = count+1;
            end
            
          end
          fprintf(['Done!\n']);
          
          obj.b = bs;
          obj.W = Ws;
          
        case 'smooth'
          bs = cell(1,d.T-burn_in);
          Ws = cell(1,d.T-burn_in);
          NT = length(burn_in:p.downsample_est:d.T);
          bs0 = cell(1,NT);
          Ws0 = cell(1,NT);
          obj.CIF = zeros(d.T - burn_in, 1);
          X_bwd = flipud(obj.X);
          y_bwd = fliplr(obj.y);
          [b_bwd,stats_bwd] = obj.glmfit0(X_bwd,y_bwd,obj.link);
          W_bwd = stats_bwd.covb;
          
          % noise (Kalman gain) matrix:
          noise_mtx = zeros(N_cov);
          for n = 1:N_cov_types
            for i = p.covariatee_ind{n}
              noise_mtx(i, i) = p.noise(n);
            end
          end
          noise_mtx = noise_mtx(fs_update_ind,fs_update_ind);
          
          fprintf(['Filtering backward...']);
          for t = 1:d.T-burn_in
            dL_dB = X_bwd(t,fs_update_ind);
            lt = exp(X_bwd(t,:)*b_bwd);
            A = W_bwd(fs_update_ind,fs_update_ind)+noise_mtx;
            U = dL_dB';
            V = dL_dB;
            C0 = lt;
            C1 = (1/C0+V*A*U);
            % matrix inversion lemma
            W_bwd(fs_update_ind,fs_update_ind) = A-A*U*(C1\V)*A;
            dB = W_bwd(fs_update_ind,fs_update_ind)*dL_dB'*(obj.y(t) - lt);
            b_bwd(fs_update_ind) = b_bwd(fs_update_ind) + dB;
          end
          b = b_bwd; W = W_bwd;
          fprintf(['Done!\n']);
          
          % or, instead of filtering backward,
          % initialize using glmfit
%           [b,dev,stats] = glmfit(obj.X,obj.y,'poisson','constant','off');
%           W = stats.covb;
          
          fprintf(['Filtering forward...']);
          for t = 1:d.T-burn_in
            dL_dB = obj.X(t,fs_update_ind);
            lt = exp(obj.X(t,:)*b);
            A = W(fs_update_ind,fs_update_ind)+noise_mtx;
            U = dL_dB';
            V = dL_dB;
            C0 = lt;
            C1 = (1/C0+V*A*U);
            % matrix inversion lemma:
            W(fs_update_ind,fs_update_ind) = A-A*U*(C1\V)*A;
            dB = W(fs_update_ind,fs_update_ind)*dL_dB'*(obj.y(t) - lt);
            b(fs_update_ind) = b(fs_update_ind) + dB;
            obj.CIF(t) = exp(obj.X(t,:)*b);            
            bs{t} = b;
            Ws{t} = W;
          end                              
          fprintf(['Done!\n']);
          
          fprintf(['Smoothing...']);
          b = bs{end}; W = Ws{end};
          bs0{end} = b; Ws0{end} = W;
          count = NT - 1;
          for t = d.T-burn_in-1:-1:1
            W0 = Ws{t}(fs_update_ind,fs_update_ind);
            W1 = W0+noise_mtx;
            iW1 = inv(W1);
            b(fs_update_ind) = bs{t}(fs_update_ind) + W0*iW1*(b(fs_update_ind) - bs{t}(fs_update_ind));
            W(fs_update_ind,fs_update_ind) = W0 + W0*iW1*(W(fs_update_ind,fs_update_ind) - W1)*iW1*W0;
            lt = exp(obj.X(t,:)*b);
            
            if mod(t, p.downsample_est)==0
              bs0{count} = b;
              Ws0{count} = W;
              count = count-1;
            end
          end
          fprintf(['Done!\n']);
          
          obj.b = bs0;
          obj.W = Ws0;
      end
      fprintf(['Done!\n']);

      obj = obj.calcGOF(); % goodness-of-fit
    end
    
    function obj = calcGOF(obj, p)      
      switch obj.link
        case 'log'
          obj.LL = sum(log(poisspdf(obj.y,obj.CIF)));
          obj.dev = 2*(sum(log(poisspdf(obj.y,obj.y))) - obj.LL);
        case 'logit'
          obj.LL = sum(log(binopdf(obj.y,ones(size(obj.y)),obj.CIF)));
          obj.dev = 2*(sum(log(binopdf(obj.y,ones(size(obj.y)),obj.y))) - obj.LL);
      end
      obj.AIC = obj.dev+2*size(obj.b,1);

      % time rescaling & KS test (NOTE: Should only be for poisson
      % regression)
      spike_ind = find(obj.y);         
      numISIs = length(spike_ind)-1;
        
        if numISIs<3, error('Too few data points for goodness-of-fit'); end
        
        z = zeros(1, numISIs);
        for j=1:numISIs                                                
          z(j) = sum(obj.CIF(spike_ind(j)+1:spike_ind(j+1)));
        end
        
        switch obj.rs
          case 'idenity'
            rs_cdf = @(z)(expcdf(z));
          case 'exp'
            z = 1-exp(-z);
            rs_cdf = @(z)(unifcdf(z,0,1));
        end
        
        [eCDF,xCDF] = ecdf(sort(z));
        aCDF = rs_cdf(xCDF);
        ks_stat = max(abs(aCDF-eCDF));
        ks_ci = 1.96/sqrt(numISIs+1);
      
%       catch
%           z = [];
%           ks_stat = NaN;
%           ks_ci = NaN;
%       end
      
      obj.rsISI = z;
      
      % compute p-value for ks-statistic (from kstest2.m):
      % WARNING: still under construction
      n1     =  length(aCDF);
      n2     =  length(eCDF);
      n      =  n1 * n2 /(n1 + n2);
      lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * ks_stat , 0);
% % %       
% % %       %
% % %       % 1-sided test:
% % %       % % %    ks_p  =  exp(-2 * lambda * lambda);
% % %       % 2-sided test (default):
% % %       %  "Use the asymptotic Q-function to approximate the 2-sided P-value."
      j       =  (1:101)';
      ks_p  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
      ks_p  =  min(max(ks_p, 0), 1);
% % %       ks_p = NaN;
      obj.KS = [ks_stat,ks_ci,ks_p];
      
    end
    
    function obj = makeX(obj,d,p)
      for i = 1:length(p.covariate_names)
        channels = p.covariate_channels{i};
        basis = p.covariate_bases{i};     
        knots = p.covariate_knots{i};
        ind = p.covariate_ind{i};        
        Xi = obj.makeX_block(d, channels, basis, knots, p.s);
        obj.X(:,ind) = Xi; clear Xi;
      end
%       obj.X = obj.X(p.get_burn_in()+1:end,:);
    end
    
    function Xi = makeX_block(obj, data, channels, basis, knots, s)
      if nargin<6, s = 0.5; end
      s_coeff = [-s  2-s s-2  s; 2*s s-3 3-2*s -s; ...
           -s   0   s   0;   0   1   0   0];
      N = length(knots);
      
      % make X for firing rate covariate
      if isequal(channels,0)
        bins_per_knot = round(diff(knots)*data.T);
        count = 0;
        
        switch basis
          case 'spline'
            Xi = zeros(data.T, N+2);
            for n = 1:N-1
              temp=bins_per_knot(n);
              alphas=1/temp*(1:temp);
              Xi(count + (1:temp), n+(0:3)) = Xi(count + (1:temp), n+(0:3)) + ...
                  [alphas'.^3 alphas'.^2 alphas' ones(temp, 1)]*s_coeff;
              count=count+bins_per_knot(n);
            end
          case 'indicator'
            Xi = zeros(data.T, N-1);
            for n = 1:N-1
              Xi(count+(1:bins_per_knot(n)),n) = 1;
              count = count+bins_per_knot(n);
            end
        end
        
      % make X for other covariates
      else
        
        if length(channels)<=1
          d = data.dn(channels,:);
        else
          d = sum(data.dn(channels,:),1);
        end 
        
        switch basis
          case 'spline'
            if range(knots)>1, dtau = 1; else dtau = 0.01; end;
            tau = knots(1):dtau:knots(end);
            NT = length(tau);
            onset = knots(1);
            offset = knots(end);
            intvls = diff(knots)*1/dtau;                        
            
            Xi = zeros(data.T, N+2);
            Xs = zeros(NT,N+2);
            
            count=1;
            for n=1:N-1
              I = intvls(n); % length of interval (num. of bins)
              alphas = (0:I-1)./I;
              Xs(count:count+I-1, n+(0:3)) = [alphas'.^3 alphas'.^2 alphas' ones(I,1)] * s_coeff;
              count = count+I;
            end
            Xs(end, N-1:N+2) = [1 1 1 1] * s_coeff; % alpha = 1
            
            spks = find(d);
            for i = 1:length(spks)
                bins_to_end = data.T - spks(i);
                if spks(i) + offset > data.T
                    Xi(spks(i)+(onset:bins_to_end),:) = ...
                        Xi(spks(i)+(onset:+bins_to_end),:) + d(spks(i)) * Xs(1:bins_to_end-onset+1,:);
                else
                    Xi(spks(i)+(onset:offset),:) = ...
                        Xi(spks(i)+(onset:offset),:) + Xs;
                end
            end
          case 'indicator'
            Xi = zeros(data.T, N);
            burn_in = knots(end);
            for n = 1:N
              Xi(burn_in+1:end, n) =  d((burn_in+1:end)-knots(n))';
            end
        end
      end
    end
    
    function perf_pred(obj)
      % function to handle perfect predictors
      %
      %
      M = obj.b;
      (num2str(M))
      
    end
    
    function [b,stats] = glmfit0(obj, X, y_in, link)
      
      if size(y_in,2)>size(y_in,1), y = y_in';
      else y = y_in; end;

      switch link
        case 'identity'
          linkFn = @(x) x;
          ilinkFn = @(x) x;
          linkFnprime = @(x) ones(size(x));
          mu = y; % initial conditions
        case 'log'
          linkFn = @(x) log(x);
          ilinkFn = @(x) exp(x);
          linkFnprime = @(x) 1./x;     
          sqrtvarFn = @(x) sqrt(x);
          mu = y + 0.25; % initial conditions
      end
      
      N = size(X,1);
      p = size(X,2);
      pwts = ones(N,1);
      b = ones(p,1);
      eta = linkFn(mu);
      
      % convergence parameters
      eps = 1e-6;
      iterLim = 100;
      offset = 1e-3;

      for iter = 1:iterLim
        z = eta - offset + (y - mu) .* linkFnprime(mu);
        b_old = b;
        deta = linkFnprime(mu);
        sqrtirls = abs(deta) .* sqrtvarFn(mu);
        sqrtw = sqrt(pwts) ./ sqrtirls;

        %  from wfit.m:
        zw = z .* sqrtw;
        Xw = X .* sqrtw(:,ones(1,p));
        [Q,R0] = qr(Xw,0);
        if iter>20
          0;
        end
        
        if ~isempty(regexp(lastwarn(),'singular','once'))
          iter
          fprintf('\nWarning: Singular matrix encountered during QR decomposition\n');
          break;
        end
        
        R = R0;
        b = R \ (Q'*zw);

        % show estimate at each step
        % [[iter;0] b]
        
        %  stop if converged:
        if norm(b - b_old, inf) < eps, break; end
        
        % stop if singular warning received:
        % DELETED
        
        eta = offset + X*b;
        mu = ilinkFn(eta);
      end
      
      % glmfit covariance:
      RI = R\eye(p);
      C = RI * RI';
      % assumes no over/under-dispersion (s=1):
      % C=C*s^2;

      stats.beta = b;
      stats.dfe = N-p;
      stats.sfit = [];
      stats.covb = C;
      stats.s = 1;
      stats.estdisp = 0;
      stats.se = sqrt(diag(stats.covb));
      stats.t = b ./ stats.se;
      stats.p = 2 * normcdf(-abs(stats.t));
      % stats.wts = diag(W);
    end

    function plot(obj, d, p)
      
      global PLOT_COLOR;
      global DO_CONF_INT;
      global DO_MASK;
      
      Z = 2;
      N_covar_types = length(p.covariate_names);
      burn_in = p.get_burn_in();
      % NOTE: ASSUMES 'rate' is first covariate
      % and all other types ('self-hist', 'ensemble', etc)
      % come afterwards
      
      % RATE----------------------
      subplot(N_covar_types,1,1); hold on;
      T0 = length(p.covariate_knots{1});
      ind = p.covariate_ind{1};
      switch obj.fit_method
        case 'glmfit'
          b1 = obj.b(ind);
          if DO_CONF_INT, W1 = obj.W(ind,ind); end;
        case {'filt','smooth'}
          b1 = obj.b{1}(ind);
          if DO_CONF_INT, W1 = obj.W{1}(ind,ind); end;
      end
      
      switch p.covariate_bases{1}
        case 'spline'
          if DO_CONF_INT
            [t_axis,Y,Ylo,Yhi] = plot_spline(p.covariate_knots{1},b1,p.s,W1,Z);
            t_axis = t_axis*(d.t(end)-d.t(1)) + d.t(1); % convert to secs
            L = exp(Y')/d.dt; Llo = exp(Ylo')/d.dt; Lhi = exp(Yhi')/d.dt;
            shadedErrorBar(t_axis,L,[Lhi-L; L-Llo],{'Color',PLOT_COLOR});
          else
            [t_axis,Y] = plot_spline(p.covariate_knots{1},obj.b(ind),p.s);
            t_axis = t_axis*(d.t(end)-d.t(1)); % convert to secs
            plot(t_axis,exp(Y')/d.dt,'Color',PLOT_COLOR,'linewidth',2);
          end
        case 'indicator'
          t_axis = p.covariate_knots{1} * (d.t(end)-d.t(1)) + d.t(1);
          if DO_CONF_INT
            Y = b1;
            % NOTE: NEED TO MODIFY HOW Ylo, Yhi
            % are computed
            Ylo = Y - Z*sqrt(diag(W1)); Yhi = Y + Z*sqrt(diag(W1)); % 1d case
            L = exp(Y')/d.dt; Llo = exp(Ylo')/d.dt; Lhi = exp(Yhi')/d.dt;
            for t = 1:T0-1
              shadedErrorBar([t_axis(t),t_axis(t+1)],L(t)*ones(1,2),[(Lhi(t)-L(t))*ones(1,2); (L(t)-Llo(t))*ones(1,2)],{'Color',PLOT_COLOR}); 
              hold on;
            end
          else
            for t = 1:T0-1
              plot([t_axis(t),t_axis(t+1)],exp(b1(t))/d.dt*ones(1,2),PLOT_COLOR,'linewidth',2);
            end
          end
      end
      xlabel('time [s]');
      ylabel('[Hz]');
      title(p.covariate_names{1});
      
      % OTHER COVARIATES----------
      for covar_num = 2:N_covar_types        
        subplot(N_covar_types,1,covar_num); hold on;
        ind = p.covariate_ind{covar_num};
        switch obj.fit_method
          case 'glmfit'
            switch p.covariate_bases{covar_num}
              case 'spline'
                if DO_CONF_INT
                  [lag_axis,Y,Ylo,Yhi] = plot_spline(p.covariate_knots{covar_num},obj.b(ind),p.s,obj.W(ind,ind),Z);
                  lag_axis = lag_axis*d.dt*1e3; % convert from bins to ms
                  L = exp(Y'); Llo = exp(Ylo'); Lhi = exp(Yhi');
                  shadedErrorBar(lag_axis,L,[Lhi-L; L-Llo],{'Color',PLOT_COLOR});
                else
                  [lag_axis,Y] = plot_spline(p.covariate_knots{covar_num},obj.b(ind),p.s);
                  lag_axis = lag_axis*d.dt*1e3; % convert from bins to ms
                  plot(lag_axis,exp(Y'),PLOT_COLOR,'linewidth',2);
                end
              case 'indicator'
                lag_axis = p.covariate_knots{covar_num};
                if DO_CONF_INT                  
                  Y = exp(obj.b(ind));                  
                  error('write more code!');
                  % assign Y, Ylo, Yhi based on covariance structure
                  L = exp(Y'); Llo = exp(Ylo'); Lhi = exp(Yhi');
                  shadedErrorBar(lag_axis,L,[Lhi-L; L-Llo],{'Color',PLOT_COLOR});
                else
                  plot(lag_axis,exp(obj.b(ind)'),PLOT_COLOR,'linewidth',2);                  
                end
            end
            plot([lag_axis(1) lag_axis(end)],[1 1],'k--');
            xlabel('lag time [ms]');
            ylabel('mod.');

          case {'filt', 'smooth'}
            t_ind = burn_in+1:p.downsample_est:d.T; % modified time axis
            NT = length(t_ind);
            switch p.covariate_bases{covar_num}
              case 'spline'                
                lag_axis = p.covariate_knots{covar_num}(1):p.covariate_knots{covar_num}(end);
                lag_axis = lag_axis*d.dt*1e3; % covert from bins to ms
                all_L = zeros(length(lag_axis), NT);            
                for t = 1:NT
                  if DO_MASK
                    [~,Y,Ylo,Yhi] = plot_spline(p.covariate_knots{covar_num},obj.b{t}(ind),p.s,obj.W{t}(ind,ind),Z);
                    % mask y:
                    good_ind = [find(Yhi>0)', find(Ylo<0)'];
                    bad_ind = setdiff(1:length(lag_axis), good_ind);
                    Y(bad_ind) = 0;
                    all_L(:,t) = exp(Y');
                  else
                    [~,Y] = plot_spline(p.covariate_knots{covar_num},obj.b{t}(ind),p.s);
                    all_L(:,t) = exp(Y');                
                  end
                end
              case 'indicator'
                lag_axis = p.covariate_knots{covar_num}(1):p.covariate_knots{covar_num}(end);
                lag_axis = lag_axis*d.dt*1e3; % covert from bins to ms
                all_L = zeros(length(lag_axis), NT);
                for t = 1:NT
                  if DO_MASK
                    Y = exp(obj.b{t_ind(t)}(ind));
                    % assign Ylo, Yhi
                    error('write more code!');                    
                    % mask Y:
                    good_ind = [find(Yhi>0), find(Ylo<0)];
                    bad_ind = setdiff(1:length(lag_axis), good_ind);
                    Y(bad_ind) = 0;
                    all_L(:,t) = exp(Y');
                  else
                    all_L(:,t) = exp(obj.b{t_ind(t)}(ind));
                  end
                end
            end
            imagesc(d.t(t_ind), lag_axis, all_L);
            xlabel('time [s]');
            ylabel('lag time [ms]');            
        end
        xlim(round([p.covariate_knots{covar_num}(1),p.covariate_knots{covar_num}(end)*0.8]*d.dt*1e3));
        title(p.covariate_names{covar_num});
      end
      
      update_fig(); % change font size, interpreter
    end

    function gof_plot(obj, d)      
%
%	pp_model.gof_plot(d)
%
%	pp_gof.m
%	Part of the Point Process Toolbox
%	By Grant Fiddyment, Boston University, September 2012
%	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	INPUTS:-------------------------------
%	d: a point process data object
%
%	DESCRIPTION:---------------------------
%   Visualizes the goodness of fit of a point process model based on
%   several features of the conditional intensity function(/process):
%   empirical distribution of rescaled ISIs, KS-plot comparing CDF of
%   rescaled ISIs against theoretical Exp[1] CDF, QQ-plot comparing
%   empirical and theoretical distributions by their quantiles, 
%   recorded from epileptic patients at MGH. Uses the time rescaling
%   theorem to compare the rescaled ISIs to an Exp[1] distribution
%   as well as calculates their (partial/) autocorrelation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
 
      %   Conditional intensity function
      subplot(2,3,1);
      if nargin<2, obj.cif_plot(); else obj.cif_plot(d.t); end

      %   Residual process
      subplot(2,3,2);
      if nargin<2, obj.res_plot(); else obj.res_plot(d.t); end
      
      %   Rescaled ISI dist.
      subplot(2,3,3); obj.isi_plot();

      %   KS plot (Theoretical vs actual dist. of rescaled ISI)
      subplot(2,3,4); obj.ks_plot();
      
      %   QQ plot
      subplot(2,3,5); obj.qq_plot();

      %   Autocorrelation
      subplot(2,3,6); obj.ac_plot();
     
    end
    
    function cif_plot(obj, t)
      global PLOT_COLOR
      if nargin<2
        t = 1:length(obj.CIF); 
        dt = 1;
      else
        t = t(end-length(obj.CIF)+1:end);
        dt = t(2) - t(1);
      end
      
      plot(t,obj.CIF/dt,'color', PLOT_COLOR); hold on;
      title('$\lambda_t$','interpreter','latex');
      xlabel('time [s]');
      ylabel('[Hz]');
      update_fig();
    end
      
    function [ks_stat, ks_ci] = ks_plot(obj)
      global PLOT_COLOR     
      Z = 1.96; 
      % Z = 1.63;
      % Z = 1.36;
      %   calculate ks statistic, confidence bounds
      numISIs = length(obj.rsISI);      
      if numISIs>2
        [eCDF,xCDF] = ecdf(sort(obj.rsISI));
        switch obj.rs
        case 'identity'
          mycdf = @(x)(expcdf(x,1));
        case 'exp'
          mycdf = @(x)(unifcdf(x,0,1));
        end
        aCDF = mycdf(xCDF);                    
        ks_stat = max(abs(aCDF-eCDF));
        ks_ci = Z/sqrt(numISIs+1);
      else
        ks_stat = NaN;
        ks_ci = NaN;
        fprintf('Error: Too few events in data for GoF analysis\n');
        return;
      end

      % plotting
      plot(eCDF, aCDF, 'color', PLOT_COLOR, 'LineWidth', 2); hold on;
      plot([0:0.2:1], ks_ci+[0:0.2:1], 'r-',  [0:0.2:1], -ks_ci+[0:0.2:1], ...
          'r-', [0:0.2:1], [0:0.2:1], 'r--', 'LineWidth', 3);
      set(gca,'XTick',[0:0.2:1]);
      set(gca,'YTick',[0:0.2:1]);
      title('KS plot');
      xlabel('Empirical CDF');
      ylabel('Theoretical CDF');
%       text(0.05, 0.75, ['KS stat: ', num2str(ks_stat,3)]);
%       text(0.05, 0.67, ['95%  CI: ', num2str(ks_ci,3)]); 
      update_fig();
    end
      
  
    function [qx, qz] = qq_plot(obj)
      global PLOT_COLOR
      z = sort(obj.rsISI);
      if length(z)<2
        return;
      else
        I = 0:.05:1;
        qz = quantile(z,I);
        switch obj.rs
        case 'identity'
%           qx = quantile(expcdf(I,1),I);
          qx = expinv(I,1);
        case 'exp'          
%           qx = quantile(unifcdf(I,0,1),I);
          qx = unifinv(I,0,1);
        end
        plot(qx,qz,'color',PLOT_COLOR); hold on;
        plot(0:0.2:1,0:0.2:1,'r--','LineWidth',3);
        title('QQ plot');
        xlabel('Empirical quant.');
        ylabel('Theoretical quant.');
        update_fig();
      end
    end
    
    function isi_plot(obj)
      [yh,xh]=hist(obj.rsISI);
      bar(xh,yh./length(obj.rsISI));
      title('rescaled ISIs','interpreter','latex');
      xlabel('$z_j$','interpreter','latex');
      ylabel('PDF','interpreter','latex');
      update_fig();
    end
    
    function res = res_plot(obj, t)
      global PLOT_COLOR
      sum_CIF=cumsum(obj.CIF);
      num_spks=cumsum(obj.y);
      res = num_spks - sum_CIF;
      if nargin<2, t = (1:length(res));
      else t = t(end-length(res)+1:end); end
      plot(t, res, 'color', PLOT_COLOR, 'LineWidth', 2);
      title('residual process');
      ylabel('observed - estimated');
      xlabel('time [s]');
      update_fig();
    end
    
    function ac = ac_plot(obj)
      global PLOT_COLOR
      numISIs = length(obj.rsISI);
      numLags=min(200,numISIs-1);
      [ac,lags,bounds]=autocorr2(obj.rsISI,numLags,round(0.3*numLags),2);
      plot(lags, ac, [PLOT_COLOR 'x']); hold on
      plot(lags,bounds(1),'r--','LineWidth', 2);
      plot(lags,bounds(2),'r--','LineWidth', 2);
      title('autocorrelation');
      xlabel('lags');
      update_fig();
    end
  end
end
       