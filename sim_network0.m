% sim_network.m
% Given a list of model coefficients (bs)
% and *single* set of model parameters p (i.e. assumed
% to be the same across models)
%

function sim_data = sim_network0(bs, d, p)
  
  fprintf('\n\nSimulating network of point process data...\n');

  ens_ind = getnameidx(p.covariate_names,'ensemble');
  temp = isequal(p.covariate_bases{ens_ind},'spline');
  ens0_ind = p.covariate_ind{ens_ind}(1 + temp);
  
  iterLim = 50;
  convThresh = 1;
  burn_in = p.get_burn_in();
  N_covar = p.covariate_ind{end}(end);
  N_covar_types = length(p.covariate_names);
  Xs = cell(1, N_covar_types);
  for n = 1:N_covar_types
      Xs{n} = p.splineX(n);
  end
  
  dn_iter = zeros(d.N_channels, iterLim);
  prev_dn = zeros(d.N_channels,1);
  curr_dn = zeros(d.N_channels,1);
  sim_dn = zeros(d.N_channels, d.T);
  cif = zeros(d.N_channels, d.T);
  Xit = zeros(d.N_channels, N_covar);
  
  % simulate burn-in period
  % each channel has its own rate
  for i = 1:d.N_channels
    cif(i, 1:burn_in) = exp(bs{i}(1));
    sim_dn(i, 1:burn_in) = poissrnd(cif(i,1),[1,burn_in]);
  end
  
  
  self_hist_on = p.covariate_knots{2}(1);
  self_hist_off = p.covariate_knots{2}(end);
  ens_on = p.covariate_knots{ens_ind}(1);
  ens_off = p.covariate_knots{ens_ind}(end);
%   covar_on = p.covariate_knots{n}(1);
%   covar_off = p.covariate_knots{n}(end);
  
  % time loop:
  for t = burn_in+1:d.T
    
    if mod(t,10)==0, fprintf([num2str(t) '\n']); end
    
    % compute each channel's probability
    % WITH causal effects
    % WITHOUT 0-lag ensemble effects
    for i = 1:d.N_channels
      dn_i = sim_dn(i,t-self_hist_on:-1:t-self_hist_off);
      dn_ens = sum(sim_dn([1:i-1,i+1:end],t-ens_on:-1:t-ens_off),1);
      Xit(i,:) = [1 dn_i*Xs{2} dn_ens*Xs{3}];
      % more general
      %for n = 1:N_covar_types
      %  dnn = 
      %  Xit(i,p.covariate_ind{n}) = [dnn * Xs{n}];
      %end
    end
    
    % sampling loop:
    for iter = 1:iterLim
    
      % MAYBE PERMUTE CHANNEL ORDER EACH TIME?
      
      % channel loop:      
      for i = 1:d.N_channels
        Xit(i, ens0_ind) = sum(curr_dn) - curr_dn(i);
        lambda(i) = exp(Xit(i,:)*bs{i}); % compute CIF
        curr_dn(i) = poissrnd(lambda(i)); % draw sample
      end
      
      % convergence check:
      % how much has total # of spikes / 
      % channel-wise pattern of spikes changed
      % since last iter?
      if sum(abs(prev_dn - curr_dn)) < convThresh
        break;
      end
      
      % save pattern for this iteration:
      dn_iter(:,iter) = curr_dn;      
    end
    
    % show iteration process
    % does pattern converge?
%     d2 = pp_data(dn_iter,1:iterLim); 
%     d2.plot('raster'); pause; hold off;
    
    % save final pattern for time t
    cif(:,t) = lambda;
    sim_dn(:,t) = curr_dn;
    
  end
  
  % return pp_data object
  sim_data = data(sim_dn, d.t, 'Simulated data');
  fprintf('Done!\n\n');
end
