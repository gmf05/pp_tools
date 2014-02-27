% sim_network.m
% Given a list of model coefficients (bs)
% and *single* set of model parameters p (i.e. assumed
% to be the same across models)
%

function [sim_data, cif] = sim_network(bs, d, p)
  
  fprintf('\n\nSimulating network of point process data...\n');
  burn_in = p.get_burn_in();
  L = zeros(d.N_channels, 1);
  sim_dn = zeros(d.N_channels, d.T);
  cif = zeros(d.N_channels, d.T);
 
      
  N_covar_types = length(p.covariate_names);
  Xs = cell(1, N_covar_types);
  for n = 1:N_covar_types
      Xs{n} = p.splineX(n);
  end

  self_hist_on = p.covariate_knots{2}(1);
  self_hist_off = p.covariate_knots{2}(end);
  ens_on = p.covariate_knots{3}(1);
  ens_off = p.covariate_knots{3}(end);
  
  % simulate burn-in period
  % each channel has its own rate
  for i = 1:d.N_channels
    cif(i, 1:burn_in) = exp(bs{i}(1));
    sim_dn(i, 1:burn_in) = poissrnd(cif(i,1),[1,burn_in]);
  end

  % time loop:
  for t = burn_in+1:d.T
    
    if mod(t,1000)==0, fprintf([num2str(t) '\n']); end
    
    for i = 1:5
      dn_i = sim_dn(i,t-self_hist_on:-1:t-self_hist_off);
      dn_ens = sum(sim_dn([1:i-1,i+1:end],t-ens_on:-1:t-ens_off),1);
%       dn_ens = sum(d.dn([1:i-1,i+1:end],t-ens_on:-1:t-ens_off),1)/5;
      Xit(i,:) = [1 dn_i*Xs{2} dn_ens*Xs{3}];
      L(i) = exp(Xit(i,:)*bs{i});
    end
    
    cif(:,t) = L;
    sim_dn(:,t) = poissrnd(L); % draw sample
    if any(sim_dn(:,t)>1), error('Too much spiking'); end
  end
  
  % return pp_data object
  sim_data = pp_data(sim_dn, d.t, 'Simulated data');
  fprintf('Done!\n\n');
end

