function [dn,cif] = pp_sim(m, d, p, shFlag)
% initialize variables
m.X = []; m = m.makeX(d,p);
Ncov = length(p.covariate_names);
dn = zeros(1,d.T);
cif = zeros(1,d.T);
Xt = zeros(1,p.covariate_ind{end}(end));
burn_in = p.get_burn_in();

% shFlag describes whether there are self-history components
% if not, simulation doesn't require stepping through time
% and can be vectorized for faster performance:
% (NOTE: will probably want to automatically detect/set shFlag...)

% shFlag = true;

if Ncov>1 && shFlag
  shcov = 2; % self history
  etccov = setdiff(1:Ncov,shcov); % "et cetera" covariates
  shind = p.covariate_ind{shcov}; % "self-history" indices
  etcind = [p.covariate_ind{etccov}]; % et cetera indices -- not self history
  lag1 = p.covariate_knots{shcov}(1);
  lag2 = p.covariate_knots{shcov}(end);
  X0 = p.splineX(shcov);

  % dn(1:burn_in) = poissrnd(exp(m.b(1)),[1,burn_in]);
  for t=burn_in+1:d.T
    if mod(t,1e4)==0, t, end
    % compute design matrix X at time t:
    % get previous spike time lags
    dlag = dn(t - (lag1:lag2));
    Xt(shind) = sum(X0(dlag>0,:),1); % ,1 option is VERY IMPORTANT!!!
    Xt(etcind) = m.X(t,etcind);
    % compute cif
  %   L1(t) = exp(Xt(shind)*m.b(shind));
  %   L2(t) = exp(Xt(etcind)*m.b(etcind));
    Lt = exp(Xt*m.b);
    cif(t) = Lt;
    dn(t) = (poissrnd(Lt)>0);

  %   Nspks=Nspks+dn(t);
    if Lt>1e4, error('CIF blew up'); end
  end
  
else
  disp('No self history');
  cif = exp(m.X*m.b)';
  dn = (poissrnd(cif)>0);
end

end